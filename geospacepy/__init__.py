"""
Geospacepy init file
"""
import yaml,pkg_resources
config_file_as_str = pkg_resources.resource_string(__name__,'geospacepy_config.json')

#Use yaml to import json because json module imports unicode by default which is annoying
config = yaml.safe_load(config_file_as_str)
