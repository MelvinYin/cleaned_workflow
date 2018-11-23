from utils import move_replace
import os
import shutil

folder1 = "test1"
folder2 = "test2"
os.mkdir(folder1)
os.mkdir(folder2)
move_replace()



shutil.rmtree("test1")
