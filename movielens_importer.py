# Basic imports
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib

# Load data from text files.
# Simple transformation tasks
items = pd.io.parsers.read_table("u.item", sep="|", header=None) \
    .rename(columns={0:"movieId",1:"title",2:"date",4:"url"})

data = pd.io.parsers.read_table("u.data", sep="\t", header=None) \
    .rename(columns={0:"userId",1:"movieId",2:"rating",3:"timestamp"}) \
    .set_index("timestamp") \
    .sort_index()

ds = pd.merge(data, items[["movieId", "title"]], on="movieId")
