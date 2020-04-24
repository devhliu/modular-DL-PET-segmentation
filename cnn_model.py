"""
Created on 12/16/2018

@author: Kevin H. Leung

Deep Convolutional Neural Network for tumor segmentation in PET images
Keras implementation
"""

import tensorflow as tf
import keras
import numpy as np
import h5py
import random
import scipy.io as sio
from matplotlib import pyplot as plt
from random import randint

from keras.models import Model
from keras.layers import Input, Conv2D, Conv2DTranspose, MaxPooling2D, Dropout, Add
from keras.layers.advanced_activations import LeakyReLU
from keras.initializers import Constant
from keras.models import load_model

from keras import backend as K

K.set_image_data_format('channels_last')

#  Define loss function
def loss_fn(y_true, y_pred):
    # Adjust for class imbalances by giving more weighting to the tumor class when computing the cost
    y_true = tf.to_int32(tf.reshape(y_true, shape=[-1, 128, 128]))  # labels, shape=[batch_size, height, width]
    class_weights = tf.add(tf.to_float(tf.multiply(y_true, 2)), tf.to_float(tf.equal(y_true, 0)))
    # weighting on foreground class and on background class
    return tf.reduce_mean(tf.multiply(tf.nn.sparse_softmax_cross_entropy_with_logits(labels=y_true, logits=y_pred), class_weights))

# Define DICE as accuracy metric
def dice_coef(y_true, y_pred):
    y_pred = tf.to_float(K.argmax(y_pred, axis=3))
    y_true = tf.reshape(y_true, shape=[-1, 128, 128])
    intersection = K.sum(y_true * y_pred, axis=[1, 2])
    union = K.sum(y_true, axis=[1, 2]) + K.sum(y_pred, axis=[1, 2])
    return K.mean((2. * intersection) / union, axis=0)  # average dice over examples

# Repeating layers throughout the network 
def add_common_layers(filters, layer, bias_ct=0.03, leaky_alpha=0.01, drop_prob=0.1):
    layer = Conv2D(filters, (3, 3), # num. of filters and kernel size 
                   strides=1,
                   padding='same',
                   use_bias=True,
                   kernel_initializer='glorot_normal', 
                   bias_initializer=Constant(value=bias_ct))(layer)
    layer = LeakyReLU(alpha=leaky_alpha)(layer) # activation function
    layer = Dropout(drop_prob)(layer) 
    return layer


def get_cnn(num_classes=2):

    # num_classes - Tumor segmentation has two classes - "lesion" and "no lesion" pixels

    # This model has skip connections in place, here we use element-wise addition.

    # Define Convolutional Neural Network

    # Input shape 
    input = Input(shape=(128,128,1))

    # Conv1
    x = add_common_layers(16, input)

    # Conv2
    conv2 = add_common_layers(16, x)
    x = MaxPooling2D(pool_size=(2,2), strides=(2,2), padding='same')(conv2)

    # Conv3
    x = add_common_layers(32, x)

    # Conv4
    conv4 = add_common_layers(32, x)
    x = MaxPooling2D(pool_size=(2,2), strides=(2,2), padding='same')(conv4)

    # Conv5
    x = add_common_layers(64, x)

    # Conv6
    x = add_common_layers(64, x)

    # Conv7
    conv7 = add_common_layers(64, x)
    x = MaxPooling2D(pool_size=(2,2), strides=(2,2), padding='same')(conv7)

    # Transposed Convolution (upsampling)
    x = Conv2DTranspose(64, (2,2), strides=(2,2), padding='same', use_bias=True, kernel_initializer='glorot_normal', bias_initializer=Constant(value=0.03))(x)
    x = LeakyReLU(alpha=0.01)(x)

    # Conv8
    x = Add()([x, conv7])
    x = add_common_layers(64, x)

    # Conv9
    x = add_common_layers(64, x)

    # Conv10
    x = add_common_layers(64, x)

    # Transposed Convolution (upsampling)
    x = Conv2DTranspose(32, (2,2), strides=(2,2), padding='same', use_bias=True, kernel_initializer='glorot_normal', bias_initializer=Constant(value=0.03))(x)
    x = LeakyReLU(alpha=0.01)(x)

    # Conv11
    x = Add()([x, conv4])
    x = add_common_layers(32, x)

    # Conv12
    x = add_common_layers(32, x)

    # Transposed Convolution (upsampling)
    x = Conv2DTranspose(16, (2,2), strides=(2,2), padding='same', use_bias=True, kernel_initializer='glorot_normal', bias_initializer=Constant(value=0.03))(x)
    x = LeakyReLU(alpha=0.01)(x)

    # Conv13
    x = Add()([x, conv2])
    x = add_common_layers(16, x)

    # Conv14
    x = Conv2D(num_classes, (3,3), strides=1, padding='same', use_bias=True, kernel_initializer='glorot_normal', bias_initializer=Constant(value=0.03))(x)
    x = LeakyReLU(alpha=0.01)(x)

    output = x

    model = Model(inputs=[input], outputs=[output])

    model.compile(loss=loss_fn, optimizer='adam', metrics=[dice_coef]) 

    return model
    
