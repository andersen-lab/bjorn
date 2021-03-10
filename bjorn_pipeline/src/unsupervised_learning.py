from PIL import Image
import matplotlib.pylab as plt
import numpy as np
from scipy.spatial import distance
from scipy.stats import mode
from sklearn.datasets import load_breast_cancer
from sklearn.ensemble import RandomForestClassifier
from sklearn.cluster import SpectralClustering
from sklearn.metrics import confusion_matrix, accuracy_score
from skimage.io import imread


def spectral_clustering_unsupervised(X, num_classes=2, n_trees=50):
    """Given dataset of (n_samples, m_features)
    Return predicted classifications (num_classes)
    using Brieman's Trick for Unsupervised Learning
    via Spectral Clustering"""
    N = len(X)  # number of samples
    X_new, y_new = create_breimans_dataset(X)  # pseudo-labeled data
    # Fit a random forest (RF) model on the new dataset
    rf = RandomForestClassifier(n_estimators=n_trees, oob_score=True)
    rf.fit(X_new, y_new);
    # Get sample indexes in each terminal node of the RF model
    sample_ids_per_leaf = get_leaf_samples(rf, X)
    # Compute similarity matrix between the samples 
    S = compute_affinity_matrix(sample_ids_per_leaf, N)
    # Perform Spectral Clustering
    cluster_model = SpectralClustering(n_clusters=num_classes, 
                                       affinity='precomputed')
    preds = cluster_model.fit_predict(S)
    # WARNING: predictions could be flipped! i.e. 1-0 instead of 0-1 (fix this)
    return preds



def create_breimans_dataset(X):
    """Create dataset of X and X_scrambled with 
    pseudo-labels 0 and 1, respectively"""
    N = X.shape[0]
    X_scramble = np.random.permutation(X)
    y_normal = np.zeros(N)
    y_scramble = np.ones(N)
    X_new = np.concatenate([X, X_scramble])
    y_new = np.concatenate([y_normal, y_scramble])
    return X_new, y_new

"""Tests
assert X_new.shape[0] == X.shape[0]*2
"""


def kmeans_img(img_fname:str, k=4, centroids='kmeans++', 
               tolerance=0.01, color=False):
    """Compress Image using K-Means Clustering"""
    
    if color: 
        img = imread(img_fname)
        h, w, c = img.shape
        X = np.reshape(img, (w * h, c)).astype(np.uint8)
    else:
        img = Image.open(img_fname)
        img = img.convert("L")     # grayscale
        w, h = img.size
        X = np.array(img.getdata(), dtype=np.uint8)

    centroids, clusters = kmeans(X, k=k, 
                                 centroids='kmeans++', 
                                 tolerance=0.01)
    if color: 
        X = reassign_colors(X, centroids, clusters)
        img_ = Image.fromarray(X.reshape((h,w,c)))
    else: 
        X = reassign_colors(X, centroids, clusters)
        img_ = Image.fromarray(X.reshape(h,w), 'L')
    return img_


def kmeans(X:np.ndarray, k:int, centroids=None, tolerance=1e-2):
    """K-Means Clustering"""
    # Initialize Input, Centroids, Clusters
    if len(X.shape) < 2: X = np.expand_dims(X, axis=1)
    if centroids: 
        centroids = init_centroids(X, k)
    else:
        centroids = X[np.random.choice(X.shape[0], size=k, replace=False), :]
    clusters = [[] for c in centroids]
    # Until Convergence
    while(True):
        # Centroid Norm at t-1
        old_centroids_norm = np.linalg.norm(centroids)
        # Minimum Distance to Centroid Per Record
        js = np.argmin(distance.cdist(X, centroids), axis=1)
        # Assign Records to Clusters
        for i, c in enumerate(clusters):
             clusters[i] = np.where(js == i)
        # Recompute Cluster Centroids
        for j, c in enumerate(centroids):
            centroids[j] = np.mean(X[clusters[j]], axis=0)
        # Centroid Norm at t
        centroids_norm = np.linalg.norm(centroids)
        # Convergence!
        if np.abs(centroids_norm - old_centroids_norm) < tolerance:
            preds = np.zeros_like(X)
            clusters = np.asarray(clusters)[np.argsort(centroids.squeeze())[::-1]]
            for i, c in enumerate(clusters):
                preds[c] = i
            return preds

        
def get_preds_from_clusters(clusters, N):
    """Compute (and spit out) Accuracy and Confusion Matrix
    Used for Classification tasks using K-Means Clustering"""
    
    preds = np.zeros_like(N)
    for i, c in enumerate(clusters):
        preds[c] = i
    return preds
        

def reassign_colors(X, centroids, clusters):
    """Reassigns pixel values based on Clusters and Centroids. 
    Used in Image Compression using K-Means"""
    centroids = centroids.astype(np.uint8)
    X_ = np.zeros_like(X)
    for i, c in enumerate(clusters):
        X_[c] = centroids[i]
    return X_


def init_centroids(X, k):
    """Initialize Centroids using K-Mean++"""
    centroids = np.zeros((k, X.shape[1]))
    centroid_idxs = [np.random.choice(X.shape[0], size=1)[0]]
    centroids[0] = X[centroid_idxs[0], :]
    for j in range(1, k):
        centroid_idxs = get_next_centroid(X, centroid_idxs)
        centroids[j] = X[centroid_idxs[-1], :]
    return centroids


def get_next_centroid(X, centroid_idxs):
    """Compute New Centroid as part of K-Means++
    Initialization"""
    X_idxs = np.arange(X.shape[0])
    # remove idxs in centroid_idxs 
    mask = np.isin(X_idxs, centroid_idxs)
    dists = distance.cdist(X[~mask], X[centroid_idxs, :])
    # compute min distance for each record
    min_dists = np.min(dists, axis=1)
    # pick record idx with max min distance
    max_dist_idx = np.argmax(min_dists)
    # recover original idx of that record
    max_dist_idx = X_idxs[~mask][max_dist_idx]
    # update list of initial centroids
    centroid_idxs.append(max_dist_idx)
    return centroid_idxs


def likely_confusion_matrix(targets, clusters):
    """Compute (and spit out) Accuracy and Confusion Matrix
    Used for Classification tasks using K-Means Clustering"""
    preds = np.zeros_like(targets)
    for c in clusters:
        preds[c] = mode(targets[c])[0][0]
    cm = confusion_matrix(targets, preds)
    acc = accuracy_score(targets, preds)
    print(f"Confusion matrix: \n {cm}")
    print(f"clustering accuracy: {acc:.4f}")


def get_leaf_samples(rf, X:np.ndarray):
    """Given a trained RandomForest model (rf) and data (X)
    Return the indexes of the data samples in each terminal node"""
    n_trees = len(rf.estimators_)
    leaf_samples = []
    leaf_ids = rf.apply(X)
    for t in range(n_trees):
        uniq_ids = np.unique(leaf_ids[:, t])
        sample_idxs_per_leaf = [np.where(leaf_ids[:, t] == idx)[0] \
                                for idx in uniq_ids]
        leaf_samples.extend(sample_idxs_per_leaf)
    return leaf_samples

def compute_affinity_matrix(sample_idxs_per_leaf, N):
    """Given sample indexes in each terminal node
    Compute a similarity matrix (n_samples, n_samples)"""
    S = np.zeros((N, N))
    for sample_ids in sample_idxs_per_leaf:
        if len(sample_ids) == 1: continue
        for idx in sample_ids:
            mask = np.setdiff1d(sample_ids, idx)
            S[idx][mask] += 1
    return S


def train_rf_model(X_new, y_new, n_trees=50):
    rf = RandomForestClassifier(n_estimators=n_trees,
                                oob_score=True)
    rf.fit(X_new, y_new);
    return rf


def visualize_pair(img1, img2, k=4, figsize=(15, 30)):
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].set_title(f"Compressed using {k} Clusters")
    ax[0].set_title("Original")
    if len(img1.shape) > 2:
        ax[0].imshow(img1)
        ax[1].imshow(img2)
    else:
        ax[0].imshow(img1, cmap=plt.cm.gray)
        ax[1].imshow(img2, cmap=plt.cm.gray)
    plt.show()
    

def visualize_multi(imgs, ks, figsize=(15, 30)):
    f, axes = plt.subplots(2, 4, figsize=figsize)

    axes[0][0].set_title("Original")
    axes[0][0].imshow(imgs[0])
    axes[0][1].set_title(f"{ks[0]} Clusters")
    axes[0][1].imshow(imgs[1])
    axes[0][2].set_title(f"{ks[1]} Clusters")
    axes[0][2].imshow(imgs[2])
    axes[0][3].set_title(f"{ks[2]} Clusters")
    axes[0][3].imshow(imgs[3])
    axes[1][0].set_title(f"{ks[3]} Clusters")
    axes[1][0].imshow(imgs[4])
    axes[1][1].set_title(f"{ks[4]} Clusters")
    axes[1][1].imshow(imgs[5])
    axes[1][2].set_title(f"{ks[5]} Clusters")
    axes[1][2].imshow(imgs[6])
    axes[1][3].set_title(f"{ks[6]} Clusters")
    axes[1][3].imshow(imgs[7])
    plt.show()
    

def visualize_multi2(imgs, ks, figsize=(15, 30)):
    f, axes = plt.subplots(2, 2, figsize=figsize)

    axes[0][0].set_title("Original")
    axes[0][0].imshow(imgs[0])
    axes[0][1].set_title(f"{ks[0]} Clusters")
    axes[0][1].imshow(imgs[1])
    axes[1][0].set_title(f"{ks[1]} Clusters")
    axes[1][0].imshow(imgs[2])
    axes[1][1].set_title(f"{ks[2]} Clusters")
    axes[1][1].imshow(imgs[3])
    plt.show()