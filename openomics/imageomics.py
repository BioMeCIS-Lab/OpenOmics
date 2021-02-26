import os

import h5py
import large_image
import numpy as np
from dask import delayed

# import histomicstk as htk
# import histomicstk.segmentation.positive_pixel_count as ppc


class WholeSlideImage:
    def __init__(self, cohort_name, folder_path, force_preprocess=False):
        """
        Args:
            cohort_name:
            folder_path:
            force_preprocess:
        """
        self.cancer_type = cohort_name
        if not os.path.isdir(folder_path) or not os.path.exists(folder_path):
            raise NotADirectoryError(folder_path)

        fname = os.path.join(folder_path, "models", "wsi_preprocessed.hdf5")

        f = h5py.File(fname, "w")

        if (not "wsi_preprocessed" in f) or force_preprocess:
            print("Preprocessing new WSI's")
            self.run_preprocess(f, folder_path)

        else:
            print("Already has wsi_preprocessed. Loading data from hdf5 file")

    @classmethod
    def name(cls):
        return __class__.__name__

    def run_preprocess(self, f, folder_path):
        """
        Args:
            f:
            folder_path:
        """
        wsi_preprocessed = f.create_dataset("wsi_preprocessed", (100, ),
                                            dtype="i")
        wsi_file = self.wsi_file_iterator(folder_path)

        i = 2
        while True and i > 0:
            imagePath = os.path.join(folder_path, wsi_file.__next__())
            i = i - 1
            self.preprocess_wsi(f, imagePath)

    def preprocess_wsi(self, f, imagePath):
        """
        Args:
            f:
            imagePath:
        """
        print(imagePath)
        print(slide_to_tile(imagePath))
        pass

    def wsi_file_iterator(self, folder_path):
        """
        Args:
            folder_path:
        """
        has_any_wsi = False
        for file in os.listdir(folder_path):
            if file.endswith(".svs"):
                has_any_wsi = True
                yield file

        if not has_any_wsi:
            raise Exception("Folder " + folder_path +
                            " doesn't contain any WSI .svs files")


def slide_to_tile(slide_path, params=None, region=None, tile_grouping=256):
    """Function to parallelize any function by tiling the slide. This routine
    can also create a label image.

    Args:
        slide_path (string (path)): Path to the slide to analyze.
        params (Parameters): An instance of Parameters, which see for further
            documentation
        region (dict, optional): A valid region dict (per a large_image
            TileSource.tileIterator's region argument)
        tile_grouping (int): The number of tiles to process as part of a single
            task

    Returns:
        * **stats** (*Output*) -- Various statistics on the input image. See
          Output.
        * **label_image** (*array-like, only if make_label_image is set*)

    Notes:
        The return value is either a single or a pair -- it is in either case a
        tuple. Dask is used as configured to compute the statistics, but only if
        make_label_image is reset. If make_label_image is set, everything is
        computed in a single-threaded manner.
    """
    ts = large_image.getTileSource(slide_path)
    print(ts.getMetadata())
    kwargs = dict(format=large_image.tilesource.TILE_FORMAT_NUMPY)
    if region is not None:
        kwargs["region"] = region
    else:
        results = []
        total_tiles = ts.getSingleTile(**kwargs)["iterator_range"]["position"]
        for position in range(0, total_tiles, tile_grouping):
            results.append(
                delayed(_count_tiles)(
                    slide_path,
                    params,
                    kwargs,
                    position,
                    min(tile_grouping, total_tiles - position),
                ))
        results = delayed(_combine)(results).compute()
    return results


def _count_tiles(slide_path, params, kwargs, position, count):
    """
    Args:
        slide_path:
        params:
        kwargs:
        position:
        count:
    """
    ts = large_image.getTileSource(slide_path)

    subtotal = np.array((0, 0))
    for pos in range(position, position + count):
        tile = ts.getSingleTile(tile_position=pos, **kwargs)["tile"]
        subtotal = subtotal + np.array(tile.shape[0:2])

    return subtotal


def _combine(results):
    """
    Args:
        results:
    """
    total = np.sum(results, axis=0)
    return total


if __name__ == "__main__":
    wsi = WholeSlideImage(
        "LUAD",
        "/media/jonny_admin/540GB/Research/TCGA_LUAD-WSI/",
        force_preprocess=True,
    )
