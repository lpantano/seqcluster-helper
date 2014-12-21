import os


def splitext_plus(f):
    """Split on file extensions, allowing for zipped extensions.
    copy from bcbio
    """
    base, ext = os.path.splitext(f)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext
