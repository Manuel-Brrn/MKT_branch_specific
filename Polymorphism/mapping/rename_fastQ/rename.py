import os
import sys


def split_name(name: str):
    split_ = name.split("_")
    left, *right = split_[-1].split(".")
    split_[-1] = left
    split_.append(".".join(right))
    return split_


def rename_files(path):
    for files in os.listdir(path):
        # if not os.path.isfile(files):
        #    print(f"Ignored : {files} not a file.")
        #    continue

        id_list = split_name(files)

        if len(id_list[-1]) > 0 and id_list[-1][0] == "R":
            print(f"Ignored : {files} (Correct format)")
            continue

        elif len(id_list) == 6:
            vals = "_".join((*id_list[0:3], id_list[4])) + "." + id_list[3] + "." + id_list[5]

        elif len(id_list) == 5:
            vals = "_".join(id_list[0:3]) + "." + id_list[3] + "." +  id_list[4]

        else:
            print(f"Can not proceed {files}, ({len(id_list) - 1} underscore instead of 4 or 5)")
            continue
        os.rename(f"{path}/{files}", f"{path}/{vals}")
        print(files, "-->", vals)


rename_files(sys.argv[1])
print("done")
