library(readr)
library(dplyr)
library(dtplyr)
library(tidyr)
library(stringr)
library(classInt)
library(jsonlite)
library(sf)
library(ggplot2)
library(gganimate)
library(magick)

# originally by @gkarthik; relocated from outbreak-info/biothings_covid19

args = commandArgs(trailingOnly=TRUE)
OUTPUT_DIR = args[1]

# define variables to loop over
EPI_VARS = c("confirmed_per_100k", "confirmed_rolling", "confirmed_rolling_per_100k", "confirmed_rolling_14days_ago_diff", "confirmed_rolling_14days_ago_diff_per_100k", "dead_per_100k", "dead_rolling", "dead_rolling_per_100k", "dead_rolling_14days_ago_diff", "dead_rolling_14days_ago_diff_per_100k")

# define geographic regions to loop over
GEO_CONSTANTS = tribble(
  ~id, ~epi_file,
  "countries", "countries_data.csv",
  "US_metros", "metros_data.csv",
  "US_states", "states_data.csv",
  "US_counties", "counties_data.csv"
)


# main function -----------------------------------------------------------
generateGifs = function(numColors = 9, exportGif = TRUE, returnJson = FALSE) {
  # loop over locations
  locations = GEO_CONSTANTS %>%
    rowwise() %>%
    mutate(breaks = list(processLocation(epi_file, id, numColors, exportGif)))
  location_df = locations %>% select(id, breaks) %>% unnest(cols = c(breaks))

  if(returnJson) {
    jsonlite::write_json(location_df, str_c(OUTPUT_DIR, "breaks.json"))
    json_df = toJSON(location_df)
    return(json_df)
  } else {
    return(location_df)
  }

}


# processLocation ---------------------------------------------------------
# 1. loads in the geographic shapefiles; transforms to correct projection, etc.
# 2. loads in the epidemiology data for that location
# 3. For each variable:
#     • calculates Fisher breaks for the color ramp
#     • calculates a histogram based on those breaks
#     • merges data with the geographic shape file
#     • generates and saves a .gif for each
processLocation = function(epi_file, location, numColors, exportGif = TRUE) {
  # loop over variables
  breaks = lapply(EPI_VARS, function(variable) processVariable(epi_file, location, variable, numColors, exportGif = exportGif, returnJson = FALSE))
  breaks_df = breaks %>% bind_cols() %>% mutate(location = location)
  return(breaks_df)
}

readData = function(epi_file) {
  out = tryCatch(
    {
      read_csv(str_c(OUTPUT_DIR, epi_file), col_types = cols(date = col_date(format = "%Y-%m-%d")))
    },
    error=function(cond) {
      message(paste("File does not exist:", OUTPUT_DIR, epi_file))
      message("Skipping this file \n")
      return(NA)
    },
    warning=function(cond) {
      message("File import failed with these warnings: \n")
      print(cond)
      return(NULL)
    },
    finally={
    }
  )
  return(out)
}


# processVariable ---------------------------------------------------------
# Main workhorse to calculate the breaks, histograms, and generate the gifs
processVariable = function(epi_file, location, variable, numColors, maxN = 25000, exportGif = TRUE, returnJson = TRUE, returnAll = FALSE) {
  print(str_c("processing variable ", variable, " for location ", location))

 # map = cleanMap(map_file, proj4, location)

  df = readData(epi_file)

  # data.table manipulations are faster...
  dt = lazy_dt(df) %>%
    filter(!is.na(.data[[variable]]))

  if(location != "admin0") {
    US_date_threshold = "2020-03-01"
    dt = dt %>% filter(date >= US_date_threshold)
  }

  # Classify the breaks
  #print(dt)
  domain = calcBreaks(dt, variable, numColors, maxN)
  print('!')
  if(all(!is.na(domain))) {
    break_limits = tibble(midpt = (domain + domain %>% lag())/2, lower = domain %>% lag(), upper =  domain, width = upper - lower) %>%
      filter(!is.na(midpt))

    dt = dt %>%
      mutate(fill = cut(.data[[variable]], domain))

    counts = dt %>%
      group_by(date) %>%
      do(h = calcHist(.data[[variable]], breaks = domain)) %>%
      as_tibble() %>%
      unnest(cols = c(h)) %>%
      mutate(fill = cut(midpt, domain)) %>%
      left_join(break_limits, by = "midpt")

    # geo join data. data.table faster than dplyr...
    #maps = dt %>% inner_join(map, by="location_id")  %>% as_tibble()
    # %>%
    #   filter(!is.na(date))
    # maps = maps[!is.na(date),] # remove the counties w/ no data
    #sf::st_geometry(maps) = "geometry"

    # Create the gifs
    #if(exportGif) {
    #  createGif(maps, map, domain, counts, variable, location)
    #}

    if(returnAll) {
      return(list(maps = maps, blank_map = map, breaks = domain, hist = counts))
    } else {
      if(returnJson) {
        json_breaks = tibble(!!(paste0(variable, "_breaks")) := list(domain))
        jsonlite::write_json(json_breaks, str_c(OUTPUT_DIR, "breaks-", location, "-", variable, "-", Sys.Date(), ".json"))
        return(json_breaks)
      }
      return(tibble(!!(paste0(variable, "_breaks")) := list(domain)))
    }
  }
}




# calcBreaks --------------------------------------------------------------
calcBreaks = function(df, variable, numColors, maxN, style="fisher") {
  # Maximum value to sample to calculate breaks.
  # Necessary because a classification of 280,000 elements is insanely slow.
  # from classInt: "default 3000L, the QGIS sampling threshold; over 3000, the observations presented to "fisher" and "jenks" are either a samp_prop= sample or a sample of 3000, whichever is larger"
  # Doing this manually, since this MAY exclude the min/max values, AND the larger of 10% of 280,000 is REALLY slow (I assume classInt is doing some sort of sampling + replacement), and unclear if there are benefits to getting the precise breaks

  set.seed(25)
  if(variable %in% df$vars) {
    values = df %>% pull(.data[[variable]])
    # Manual sampling of the data so things don't blow up too much.
    # making sure to add the max and min value
    if(length(values) > maxN) {
      minVal = min(values)
      maxVal = max(values)
      values = values %>% sample(maxN)
      if(! minVal %in% values) {
        values = c(values, minVal)
      }
      if(! maxVal %in% values) {
        values = c(values, maxVal)
      }
    }
    print(df)

    breaks = classIntervals(values, numColors, style=style, warnLargeN = FALSE)

    if(str_detect(variable, "_diff")) {
      # Ensure the breaks are centered at 0 if it's a difference
      midpoint = which((breaks$brks < 0 & breaks$brks %>% lead() > 0) | breaks$brks == 0)

      padLength = length(breaks$brks) - 2 * midpoint; # changes from JS code, since .js 0-indexes, while R is 1-based.
      domain = breaks$brks

      # ensure that the padding is an even number, so the limits all apply
      if(padLength %% 2) {
        padLength = padLength + 1
      }

      if(padLength < 0) {
        maxVal = max(domain)
        domain = c(domain, rep(maxVal, -1*padLength) + seq(1, by=1, length.out=(-1*padLength)))
      }
      if(padLength > 0 ) {
        minVal = min(domain)
        domain = c(rep(minVal, padLength)+ seq(1, by=1, length.out=padLength), domain)
      }
    } else {
      domain = breaks$brks
    }

    return(sort(domain))
  } else {
    print(str_c("    WARNING: variable ", variable, " is not found. Skipping calculating breaks"))
    return(NA)
  }
}



# Geoprocessing -----------------------------------------------------------
# • Projects to an appropriate projection
# • For the US, upscales Hawaii/Puerto Rico and downsizes Alaska (sorry, you're just too big) and rotates/translates to a nicer location
cleanMap = function(map_file, proj4, id) {
  # Make sure to remove empty polygons. DC disappears due to mapshaper smoothing
  # DC screws things up, since it has no polygon; filter out places without geoms
  map = sf::read_sf(map_file) %>% filter(!st_is_empty(geometry))
  # convert it to Albers equal area
  map = sf::st_transform(map, proj4)

  if(id %in% c("US_states", "US_metro", "US_counties")) {

    # Based on https://github.com/hrbrmstr/rd3albers
    # and https://r-spatial.github.io/sf/articles/sf3.html#affine-transformations
    # extract, then rotate, shrink & move alaska (and reset projection)
    rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

    if(id == "US_states") {
      alaska = map[map$location_id == "USA_US-AK",]
      AK_ctr = st_centroid(alaska$geometry)
      AK_scale = 0.5
      AK = (alaska$geometry - AK_ctr) * rot((-50*pi)/180) * AK_scale + AK_ctr + c(0500000, -5000000)

      hawaii = map[map$location_id == "USA_US-HI",]
      HI_ctr = st_centroid(alaska$geometry)
      HI_scale = 1.75
      HI = (hawaii$geometry - HI_ctr) * rot((-35*pi)/180) * HI_scale + HI_ctr + c(2.75e6, 3.5e6)

      puertorico = map[map$location_id == "USA_US-PR",]
      PR_scale = 2
      PR_ctr = st_centroid(puertorico$geometry)
      PR = (puertorico$geometry) * rot((15*pi)/180) * PR_scale + PR_ctr + c(-6.8e6,6e6)

      map = map %>% mutate(geometry = st_sfc(ifelse(location_id == "USA_US-AK", AK[1], ifelse(location_id == "USA_US-HI", HI[1], ifelse(location_id == "USA_US-PR", PR[1], geometry)))))
    }
    if(id == "US_counties") {
      # alaska <- map[map$STATEFP == "02",]
      # AK_ctr = st_centroid(alaska$geometry)
      # AK_scale = 0.5
      # AK = (alaska$geometry - AK_ctr) * rot((-50*pi)/180) * AK_scale + AK_ctr + c(0500000, -5000000)
      #
      # hawaii <- map[map$location_id == "HI",]
      # HI_ctr = st_centroid(alaska$geometry)
      # HI_scale = 1.75
      # HI = (hawaii$geometry - HI_ctr) * rot((-35*pi)/180) * HI_scale + HI_ctr + c(2.75e6, 3.5e6)
      #
      # puertorico <- map[map$location_id == "PR",]
      # PR_scale = 2
      # PR_ctr = st_centroid(puertorico$geometry)
      # PR = (puertorico$geometry) * rot((15*pi)/180) * PR_scale + PR_ctr + c(-6.8e6,6e6)
      #
      # map = map %>% mutate(geometry = st_sfc(ifelse(STATEFP == "02", AK[1], ifelse(STATEFP == "15", HI[1], ifelse(STATEFP == "72", PR[1], geometry)))))}
    }
  }
  return(map)
}


# calcHist ----------------------------------------------------------------
calcHist = function(values, breaks) {
  hist_values = hist(values, breaks = breaks, plot = FALSE)
  return(tibble(count = hist_values$counts, midpt = hist_values$mids))
}


# createGif ---------------------------------------------------------------
createGif = function(maps, blank_map, breaks, hist, variable, location) {
  fps = 4

  # Labels for histogram
  variableLabels = tibble(confirmed_per_100k = "total cases per 100,000 residents",
                          confirmed_rolling="7 day average of daily cases",
                          confirmed_rolling_per_100k = "7 day average of daily cases per 100,000 residents",
                          confirmed_rolling_14days_ago_diff = "average cases vs. 2 weeks ago",
                          confirmed_rolling_14days_ago_diff_per_100k = "average cases per 100,000 residents vs. 2 weeks ago",
                          dead_per_100k = "total deaths per 100,000 residents",
                          dead_rolling = "7 day average of daily deaths",
                          dead_rolling_per_100k = "7 day average of daily deaths per 100,00 residents",
                          dead_rolling_14days_ago_diff = "average deaths vs. 2 weeks ago",
                          dead_rolling_14days_ago_diff_per_100k = "average deaths vs. 2 weeks ago")
  geoLocations = tibble(US_states = "U.S. states", US_metros = "U.S. metropolitan areas", US_counties = "U.S. counties", admin0 = "countries")

  # Interpolate color palette
  if(str_detect(variable, "diff")) {
    colorPalette = colorRampPalette(c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
                                    space="Lab")(length(breaks) - 1)
  } else {
    colorPalette = colorRampPalette(c("#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
                                    space="Lab")(length(breaks) - 1)
  }

  # Check the histogram and map have the same number of frames
  num_frames = length(unique(hist$date))
  if(length(unique(maps$date)) != num_frames) {
    stop("Mismatch in number of frames between histogram legend and map")
  }
  # --- MAP ---
  p_map =
    ggplot(maps) +
    geom_sf(size = 0.1, aes(fill = fill, group=date)) +
    geom_sf(size = 0.2, data = blank_map, fill = NA) +
    # scale_fill_stepsn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"), limits=range(maps$breaks), breaks=maps$breaks[1:11], na.value = "white", show.limits=T, guide="colourbar") +
    scale_fill_manual(values=colorPalette, breaks = levels(maps$fill), na.value = "white", drop=FALSE) +
    labs(title = "{format(frame_time, '%d %B %Y')}") +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(size=18, hjust = 0.5)) +
    transition_time(date)

  # total line trace -----------------------------------------------------
  total = st_drop_geometry(maps) %>%
    group_by(date) %>%
    summarise(total = sum(.data[[variable]], na.rm=TRUE))

  yMax = max(total$total) * 1.1
  yMin = min(total$total) * 1.1
  xMin = min(total$date)
  xMax = max(total$date)

  p_total = ggplot(total)

  if(str_detect(variable, "diff")) {
    p_total = p_total +
      annotate(geom ="rect", xmin = xMin, xmax = xMax, ymin = 0, ymax = yMax, fill = "#fdae61", alpha = 0.25) +
      annotate(geom ="rect", xmin = xMin, xmax = xMax, ymin = 0, ymax = yMin, fill = "#abd9e9", alpha = 0.3) +
      annotate(geom="text", x = xMin, y = yMax, label = "WORSE THAN 2 WEEKS BEFORE", colour = "#f46d43", hjust = -0.025, vjust = 1.5) +
      annotate(geom="text", x = xMax, y = yMin, label = "BETTER THAN 2 WEEKS BEFORE", colour = "#4575b4", hjust = -0.025, vjust = -0.5)
  }

  p_total = p_total +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = date, y = total, group="USA"), colour = "#2c3e50", size = 1) +
    geom_point(aes(x = date, y = total, group="USA"), colour = "#2c3e50", size = 2) +
    ggtitle(str_c("Combined ", variableLabels[[variable]])) +
    scale_y_continuous(label = scales::comma) +
    theme_minimal() +
    theme(text = element_text(size=20), axis.title = element_blank(), title = element_text(size = 9)) +
    ease_aes('linear') +
    transition_reveal(date)

  # best/worst dot plot -----------------------------------------------------
  worstPlaces = st_drop_geometry(maps) %>%
    group_by(date) %>%
    mutate(rank = row_number(desc(.data[[variable]])),
           fill = cut(.data[[variable]], breaks)) %>%
    filter(rank <= 5)

  bestPlaces = st_drop_geometry(maps) %>%
    group_by(date) %>%
    mutate(rank = row_number(.data[[variable]]),
           fill = cut(.data[[variable]], breaks)) %>%
    filter(rank <= 5)

  wp = ggplot(worstPlaces) +
    geom_col(aes_string(x = "rank", y = variable, group = "name"), width=0.05, fill="#bababa") +
    geom_point(aes_string(x = "rank", y = variable, fill = "fill"), size = 2.5, shape=21) +
    geom_text(aes(x = rank, y = 0, label = name, group = name), size = 3, hjust = 1.15) +
    theme_minimal() +
    theme(
      axis.ticks.y = element_blank(),
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(size = 0.25, colour="#aabdd1"),
      title = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(15,15,1, 70), 'pt')) +
    coord_flip(clip='off') +
    scale_x_reverse() +
    scale_y_continuous(label=scales::comma, breaks = scales::pretty_breaks(n = 3)) +
    scale_fill_manual(values=colorPalette, breaks = levels(hist$fill), na.value = "white", drop=FALSE) +
    labs(title = "Worst locations") +
    transition_time(date) +
    ease_aes("linear") +
    enter_fly(y_loc = 0)

  bp = ggplot(bestPlaces) +
    geom_col(aes_string(x = "rank", y = variable, group = "name"), width=0.05, fill="#bababa") +
    geom_point(aes_string(x = "rank", y = variable, fill = "fill"), size = 2.5, shape=21) +
    geom_text(aes(x = rank, y = 0, label = name, group = name), size = 3, hjust = 1.15) +
    theme_minimal() +
    theme(
      axis.ticks.y = element_blank(),
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_line(size = 0.25, colour="#aabdd1"),
      title = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(15,15,1, 70), 'pt')) +
    coord_flip(clip='off') +
    scale_x_reverse() +
    scale_y_continuous(label=scales::comma, breaks = scales::pretty_breaks(n = 3)) +
    scale_fill_manual(values=colorPalette, breaks = levels(hist$fill), na.value = "white", drop=FALSE) +
    labs(title = "Best locations") +
    transition_time(date) +
    ease_aes("linear") +
    enter_fly(y_loc = 0)


  # --- HISTOGRAM LEGEND ---
  barWidth = min(hist %>% filter(width > 1) %>% pull(width), na.rm = TRUE) * 0.45
  maxVal = hist %>% pull(upper) %>% max()

  p_legend =
    ggplot(hist)

  if(str_detect(variable, "diff")) {
    nudge = range(hist$midpt) %>% sapply(function(x) abs(x)) %>% min() * 0.05
    p_legend = p_legend +
      geom_vline(xintercept = 0, colour = "#aabdd1", size = 0.25, linetype = 2) +
      geom_text(aes(x = 0, y = pretty(hist$count) %>% last(), label = paste("\u2190","better")), hjust = 1, nudge_x = -1*nudge, data = tibble()) +
      geom_text(aes(x = 0, y = pretty(hist$count) %>% last(), label =paste("worse", "\u2192")), hjust = 0, nudge_x = nudge, data = tibble())
  }

  p_legend = p_legend +
    geom_hline(yintercept = 0, colour = "#2c3e50") +
    geom_rect(aes(xmin=midpt - barWidth, xmax = midpt+ barWidth, ymin=0, ymax=count, fill=fill), colour = "#2c3e50", size = 0.2) +
    geom_rect(aes(ymin = -5, ymax=-2, xmin = lower, xmax = upper, fill = fill)) +
    geom_text(aes(y=-5, x=lower, label=scales::comma(round(lower/10)*10, accuracy=1)), nudge_y = -4, check_overlap = TRUE) +
    geom_text(aes(y=-5, x=maxVal %>% max(), label=scales::comma(round(maxVal/10)*10, accuracy=1)), nudge_y = -4, check_overlap = TRUE) +
    scale_fill_manual(values=colorPalette, breaks = levels(hist$fill), na.value = "white", drop=FALSE) +
    scale_y_continuous(breaks = pretty(hist$count)) +
    labs(title = paste("Number of", geoLocations[[location]]), subtitle= variableLabels[[variable]])+
    xlab(variableLabels[[variable]]) +
    enter_grow() +
    exit_shrink() +
    ease_aes('sine-in-out') +
    transition_time(date) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_text(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(size = 0.25, colour="#aabdd1"),
      axis.text.y = element_text(size = 14),
      axis.title = element_blank())

  # Create the animation frames
  map_gif = animate(p_map, fps=fps, nframes = num_frames, renderer = magick_renderer(), width = 500, height=350)
  wp_gif = animate(wp, fps=fps, nframes = num_frames, renderer = magick_renderer(), width = 150, height=125)
  bp_gif = animate(bp, fps=fps, nframes = num_frames, renderer = magick_renderer(), width = 150, height=125)
  legend_gif = animate(p_legend, fps=fps, nframes = num_frames, renderer = magick_renderer(), width = 300, height=200)
  total_gif = animate(p_total, fps=fps, nframes = num_frames, renderer = magick_renderer(), width = 500, height=126) # total height must be even

  if(length(map_gif) != length(legend_gif)) {
    stop("Mismatch in number of frames between histogram legend and map")
  }

  if(length(map_gif) != length(wp_gif)) {
    print("Mismatch in number of frames between worst places and map")
    print(length(map_gif))
    print(length(wp_gif))
    # stop("Mismatch in number of frames between worst places and map")
  }
  if(length(map_gif) != length(bp_gif)) {
    print("Mismatch in number of frames between best places and map")
    # stop("Mismatch in number of frames between best places and map")
  }

  # Combine together
  # First: zip best/worst locations
  dotplot_gif = combineGifs(wp_gif, bp_gif, num_frames)
  legend_gif_comb = combineGifs(legend_gif, dotplot_gif, num_frames, TRUE)
  combined_gif = combineGifs(legend_gif_comb, map_gif, num_frames)
  combined_gif = combineGifs(combined_gif, total_gif, num_frames, TRUE)

  # Export!
  # Note: .mp4 is ~ 200 KB while .gif is 2-4 MB so going with the smaller file.
  # image_write_gif(combined_gif, "testergif.gif", delay=1/fps)
  image_write_video(combined_gif, paste0(OUTPUT_DIR, location, "_", variable, "_", format(Sys.Date(), "%Y-%m-%d"), ".mp4"), framerate=fps)

}

combineGifs = function(gif1, gif2, num_frames, stack = FALSE, addFooter = FALSE) {
  combined_gif = image_append(c(gif1[1], gif2[1]), stack = stack)
  for(i in 2:num_frames) {
    combined = image_append(c(gif1[i], gif2[i]), stack = stack)

    if(addFooter) {
      image_read("https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/web/src/assets/logo.png") %>%
        image_scale("35") %>%
        image_extent("600x35", gravity = "northwest") %>%
        image_background("#bababa", flatten = TRUE) %>%
        # image_border("#ffffff", "17x12") %>%
        image_annotate("outbreak.info", color = "blue", size = 12,
                       location = "+10+2", gravity = "southwest")
    }

    combined_gif = c(combined_gif, combined)
  }

  return(combined_gif)
}

# invoke the function -----------------------------------------------------
# (1) Initial call: run the each variable without GIFs
# Need to loop over each variable/location
#variable = "confirmed_rolling"
#idx = 2
#breaks = processVariable(GEO_CONSTANTS$epi_file[idx], GEO_CONSTANTS$map_file[idx], GEO_CONSTANTS$proj4[idx], GEO_CONSTANTS$id[idx], variable, numColors = 9, returnJson = TRUE, exportGif = FALSE)

# (2) Final call: create GIFs, don't output the .jsons
# breaks = generateGifs()

# Alternative possibilities:
# Can also be run individually, returning a dataframe or JSON
# breaks = processVariable(GEO_CONSTANTS$epi_file[2], GEO_CONSTANTS$map_file[2], GEO_CONSTANTS$proj4[2], GEO_CONSTANTS$id[2], "confirmed_rolling", 9, returnJson = TRUE, exportGif = T)

# (Alternate 1) Can run the whole break generation in one go and save to .json
breaks = generateGifs(returnJson = TRUE, exportGif = FALSE)
