"""This colab showcases the usage of the Enformer model published in
Effective gene expression prediction from sequence by integrating long-range interactions
Žiga Avsec, Vikram Agarwal, Daniel Visentin, Joseph R. Ledsam, Agnieszka Grabska-Barwinska,
Kyle R. Taylor, Yannis Assael, John Jumper, Pushmeet Kohli, David R. Kelley"""


import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyfaidx
from kipoiseq import Interval
import kipoiseq
import gzip
import joblib
import tensorflow_hub as hub
import tensorflow.compat.v2 as tf
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Make sure the GPU is enabled
assert tf.config.list_physical_devices(
    'GPU'), 'Start the colab kernel with GPU: Runtime -> Change runtime type -> GPU'

model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file = '../../data/hg38_fasta_files/segmented_seq_chr2.fa' # replace with your fasta file.

class Enformer:

    def __init__(self, tfhub_url):
        self._model = hub.load(tfhub_url).model

    def predict_on_batch(self, inputs):
        predictions = self._model.predict_on_batch(inputs)
        return {k: v.numpy() for k, v in predictions.items()}

    @tf.function
    def contribution_input_grad(self, input_sequence,
                                target_mask, output_head='human'):
        input_sequence = input_sequence[tf.newaxis]

        target_mask_mass = tf.reduce_sum(target_mask)
        with tf.GradientTape() as tape:
            tape.watch(input_sequence)
            prediction = tf.reduce_sum(
                target_mask[tf.newaxis] *
                self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass

        input_grad = tape.gradient(prediction, input_sequence) * input_sequence
        input_grad = tf.squeeze(input_grad, axis=0)
        return tf.reduce_sum(input_grad, axis=-1)

def one_hot_encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

model = Enformer(model_path)

genome_iterator = SeqIO.parse(
    fasta_file, "fasta")

# Input: One SeqIO fasta entry
# Returns: None
# Does: Make predictions for a genomic interval. Writes predictions to specified directories.
def encode_seq(seqObj):

		# Retrieving current sample's name.
    splitTitle = seqObj.description.split(
        "chr2:134851076-136851076_")[1].split(" ")
    dirName = splitTitle[0] # Should set dirName to each sample name

    # Creating a directory to store predictions in.
    if not (os.path.exists(dirName)):
        os.system("mkdir " + dirName)
    fout = splitTitle[1]    # fout gives range of subregion

    # Predictions only provided for one 393216 bp range. This can be adjusted.
    if (fout == "chr2:135653892-136047107"):
        sequence = str(seqObj.seq).upper()
        sequence_one_hot = one_hot_encode(sequence)
        predictions = model.predict_on_batch(
        sequence_one_hot[np.newaxis])['human'][0]

        # Creating file to write predictoins to.
        filePath = dirName + "/" + fout + ".npz"
        os.system("touch " + filePath)

        # Saves prediction as a compressed numpy array (.npz).
        np.savez_compressed(filePath, predictions)

list(map(encode_seq, genome_iterator))
