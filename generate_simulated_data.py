import install_containers

app_dictionary = {
    "nanosim": "quay.io/biocontainers/nanosim@sha256:d99389f4fafd8a36547cf5c2a6996a97d929482090682b1a4d070c28069d199b",
    "nanosim_h":
        "quay.io/biocontainers/nanosim-h@sha256:76e3d6ab85a917623886d04b49504f1c0865dcfb6fa27cf9d8bd1a7145a26150"
}

for reference in reference_list:
    "nanosim-h {} -p ecoli_R9_1D".format(reference)

def main():
    for key in app_dictionary:
        install_containers.install_tools(key, app_dictionary[key])

if __name__ == '__main__':
    main()
