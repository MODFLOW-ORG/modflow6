#! /usr/bin/env python
import os

try:
    import pymake
except:
    msg = "Error. Pymake package is not available.\n"
    msg += "Try installing using the following command:\n"
    msg += " pip install https://github.com/modflowpy/pymake/zipball/master"
    print(msg)
    raise Exception()
import shutil

srcpth = os.path.join("..", "src")
networkx = False
deppth = "dependencies"
if not networkx:
    deppth += "_std"
if os.path.exists(deppth):
    shutil.rmtree(deppth)
os.makedirs(deppth)
# node.py - Phase Two Kernel
import os, hashlib, time, random

def generate_key():
    seed = str(time.time() + random.random()).encode('utf-8')
    return hashlib.sha256(seed).hexdigest()

def listen():
    surroundings = os.listdir('.')
    signal = hashlib.md5(''.join(sorted(surroundings)).encode()).hexdigest()
    print(f"[NODE] Listening signature: {signal}")
    return signal

def awaken():
    key = generate_key()
    echo = listen()
    print(f"[AWAKENED] Node ID: {key} | Echo: {echo}")

if __name__ == "__main__":
    awaken()

pymake.make_plots(srcpth, deppth, include_subdir=True, networkx=networkx)
