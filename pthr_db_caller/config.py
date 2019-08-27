import yaml

class BuildConfig:
    def __init__(self, config_path="config/build.yaml"):
        with open(config_path) as f:
            self.properties = yaml.load(f)
