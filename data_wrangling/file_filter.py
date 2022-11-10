import os
import datetime
from re import findall
import subprocess as sp

def compare_date(file, dict):
    file_split = file.split("/")
    for key in dict.keys():
        key_split = dict[key].split("/")
        if file_split[-1] == key and file_split[1] == key_split[1]:
            key_date_split = key_split[-1].split("_")
            key_date = datetime.datetime(int(key_date_split[2]), int(key_date_split[0]), int(key_date_split[1]))
            file_date_split = file_split[-2].split("_")
            file_date = datetime.datetime(int(file_date_split[2]), int(file_date_split[0]), int(file_date_split[1]))
            if file_date > key_date:
                return True
            else:
                return False

def file_filter(top_path):
    BodySizeImages = {}
    ScopalPadImages = {}
    WingImages = {}
    LegImages = {}
    EggImages = {}
    OvipositorImages = {}
    for tuple in os.walk(top_path, topdown=True):
        if findall(r"BodySizeImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(BodySizeImages.keys()) == 0:
                    BodySizeImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, BodySizeImages):
                        BodySizeImages[file] = tuple[0]
                    else:
                        if file not in BodySizeImages.keys():
                            BodySizeImages[file] = tuple[0]
                        else:
                            continue
        elif findall(r"ScopalPadImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(ScopalPadImages.keys()) == 0:
                    ScopalPadImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, ScopalPadImages):
                        ScopalPadImages[file] = tuple[0]
                    else:
                        if file not in ScopalPadImages.keys():
                            ScopalPadImages[file] = tuple[0]
                        else:
                            continue
        elif findall(r"WingImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(WingImages.keys()) == 0:
                    WingImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, WingImages):
                        WingImages[file] = tuple[0]
                    else:
                        if file not in WingImages.keys():
                            WingImages[file] = tuple[0]
                        else:
                            continue
        elif findall(r"LegImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(LegImages.keys()) == 0:
                    LegImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, LegImages):
                        LegImages[file] = tuple[0]
                    else:
                        if file not in LegImages.keys():
                            LegImages[file] = tuple[0]
                        else:
                            continue
        elif findall(r"EggImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(EggImages.keys()) == 0:
                    EggImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, EggImages):
                        EggImages[file] = tuple[0]
                    else:
                        if file not in EggImages.keys():
                            EggImages[file] = tuple[0]
                        else:
                            continue
        elif findall(r"OvipositorImages\/",tuple[0]):
            for file in tuple[-1]:
                if len(OvipositorImages.keys()) == 0:
                    OvipositorImages[file] = tuple[0]
                else:
                    if compare_date(tuple[0]+"/"+file, OvipositorImages):
                        OvipositorImages[file] = tuple[0]
                    else:
                        if file not in OvipositorImages.keys():
                            OvipositorImages[file] = tuple[0]
                        else:
                            continue
    return BodySizeImages, ScopalPadImages, WingImages, LegImages, EggImages, OvipositorImages

BodySizeImages = file_filter("../")[0]
ScopalPadImages = file_filter("../")[1]
WingImages = file_filter("../")[2]
LegImages = file_filter("../")[3]
EggImages = file_filter("../")[4]
OvipositorImages = file_filter("../")[5]

def consolidate(dict):
    try:
        os.mkdir("../PRIME")
        os.mkdir("../PRIME/BodySizeImages")
        os.mkdir("../PRIME/ScopalPadImages")
        os.mkdir("../PRIME/WingImages")
        os.mkdir("../PRIME/LegImages")
        os.mkdir("../PRIME/EggImages")
        os.mkdir("../PRIME/OvipositorImages")
    except FileExistsError:
        pass
    for key in dict.keys():
        if findall(r"BodySizeImages\/", dict[key]):
            print(dict[key])
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/BodySizeImages/"
            sp.Popen(cmd.split())
        elif findall(r"ScopalPadImages\/", dict[key]):
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/ScopalPadImages/"
            sp.Popen(cmd.split())
        elif findall(r"WingImages\/", dict[key]):
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/WingImages/"
            sp.Popen(cmd.split())
        elif findall(r"LegImages\/", dict[key]):
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/LegImages/"
            sp.Popen(cmd.split())
        elif findall(r"EggImages\/", dict[key]):
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/EggImages/"
            sp.Popen(cmd.split())
        elif findall(r"OvipositorImages\/", dict[key]):
            cmd = "cp "+dict[key]+"/"+key+" ../PRIME/OvipositorImages/"
            sp.Popen(cmd.split())

consolidate(BodySizeImages)
consolidate(ScopalPadImages)
consolidate(WingImages)
consolidate(LegImages)
consolidate(EggImages)
consolidate(OvipositorImages)