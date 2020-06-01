# modular-DL-PET-segmentation
A physics-guided modular deep-learning based automated framework for tumor segmentation in PET

# Models 
Model weights for the trained CNN in Module 2 and Module 3 are located in "models" folder. 

For example, the model weights for Module 2 can be loaded into a keras model as follows: 

from cnn_model import get_cnn
model = get_cnn()
model.load_weights('/models/module2.h5')
