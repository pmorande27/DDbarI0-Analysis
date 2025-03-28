#!/usr/bin/python
import os
def list_all(path):
    files = os.listdir(path)
    listing = []
    for file in files:
        if "Amplitude" in file:

            listing.append(file)
    return " ".join(listing)
        
