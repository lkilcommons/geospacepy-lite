from geospacepy import omnireader
import pytest
import numpy as np
from numpy import testing as nptest
import datetime,os

def test_omnireader_can_download():
	"""
	Test that we can get to the omni FTP and that
	we have a sane local directory for storing files
	"""
	dt = datetime.datetime(2005,1,1)
	od = omnireader.omni_downloader()
	od.get_cdf(dt,'5min')
	downloaded_cdf = os.path.join(od.localdir,od.filename_gen['5min'](dt))
	assert os.path.exists(downloaded_cdf)

if __name__ == '__main__':
	pytest.main()
