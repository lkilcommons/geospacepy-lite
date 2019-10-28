"""
Geospacepy init file
"""
#import yaml,pkg_resources
#config_file_as_str = pkg_resources.resource_string(__name__,'geospacepy_config')
import traceback
try:
    from geospacepy.geospacepy_config import config
except:
    print(traceback.format_exc())
    from geospacepy.default_config import config
