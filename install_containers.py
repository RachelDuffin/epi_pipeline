"""
Installs containers required in the pipeline
"""
import os
import subprocess
from config import app_dictionary


def install_tools(key, value):
    """
    Checks if docker image specified in app_dictionary already exists, if not pulls from biocontainers
    """
    if subprocess.check_output("docker images -q " + value, shell=True):
        print("INSTALLING DOCKER IMAGE: " + key + " docker image already pulled")
    else:
        print("INSTALLING DOCKER IMAGE: " + "Installing " + key)
        os.system("docker image pull " + value)
    print("----------------------------")


def main():
    for key in app_dictionary:
        value = app_dictionary[key]
        install_tools(key, value)


if __name__ == '__main__':
    main()
