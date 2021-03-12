import errno
import io
import logging
import os
from urllib.error import URLError

import dask.dataframe as dd
import requests
import sqlalchemy as sa
import validators
import astropy
from astropy.utils import data
from requests.adapters import HTTPAdapter, Retry

import openomics


@astropy.config.set_temp_cache(openomics.config["cache_dir"])
def get_pkg_data_filename(dataurl, file):
    """Downloads a remote file given the url, then caches it to the user's home
    folder.

    Args:
        dataurl: Url to the download path, excluding the file name
        file: The file path to download

    Returns:
        filename (str): A file path on the local file system corresponding to
        the data requested in data_name.
    """
    # Split data url and file name if the user provided a whole path in file_resources
    if validators.url(file):
        dataurl, file = os.path.split(file)
        dataurl = dataurl + "/"

    try:
        logging.info("Fetching file from: {}{}, saving to {}".format(
            dataurl, file, openomics.config["cache_dir"]))

        with data.conf.set_temp("dataurl", dataurl), data.conf.set_temp(
                "remote_timeout", 30):
            return data.get_pkg_data_filename(file,
                                              package="openomics.database",
                                              show_progress=True)

    except (URLError, ValueError) as e:
        raise Exception(
            "Unable to download file at {}. Please try manually downloading the files. \n{}"
            .format(os.path.join(dataurl, file), e))


def read_db(path, table, index_col):
    """
    Args:
        path:
        table:
        index_col:
    """
    engine = sa.create_engine(path)
    conn = engine.connect()
    m = sa.MetaData()
    table = sa.Table(table, m, autoload=True, autoload_with=engine)

    # conn.execute("create table testtable (uid integer Primary Key, datetime NUM)")
    # conn.execute("insert into testtable values (1, '2017-08-03 01:11:31')")
    # print(conn.execute('PRAGMA table_info(testtable)').fetchall())
    # conn.close()

    uid, dt = list(table.columns)
    q = sa.select([dt.cast(sa.types.String)]).select_from(table)

    daskDF = dd.read_sql_table(table,
                               path,
                               index_col=index_col,
                               parse_dates={"datetime": "%Y-%m-%d %H:%M:%S"})
    return daskDF


def mkdirs(outdir):
    """Make a directory.

    Args:
        outdir: directory path
    """
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass


def retry(num=5):
    """retry connection.

    define max tries num if the backoff_factor is 0.1, then sleep() will
    sleep for [0.1s, 0.2s, 0.4s, ...] between retries. It will also force a
    retry if the status code returned is 500, 502, 503 or 504.

    Args:
        num:
    """
    s = requests.Session()
    retries = Retry(total=num,
                    backoff_factor=0.1,
                    status_forcelist=[500, 502, 503, 504])
    s.mount("http://", HTTPAdapter(max_retries=retries))

    return s


def get_decompressed_text_gzip(gzip_file):
    """
    Args:
        gzip_file:
    """
    # compressedFile = StringIO()
    # compressedFile.write(gzip_file.read())
    # compressedFile.seek(0)
    return io.TextIOWrapper(gzip_file)
    # decompressedFile = gzip.GzipFile(fileobj=gzip_file, mode='rb')
    # return decompressedFile
