from driver import *
import ReedSolomon
import fieldmath
from PIL import Image
import numpy as np
from tqdm import tqdm
import math


if __name__ == '__main__':
    # init GRS to correct up to 3 errors per 30 symbols. This gives a score of 0.0107 from the driver file
    f = fieldmath.BinaryField(0x11d)
    # with [82 64 19] we get 114 it/s encode, 1.48 decode
    k = 64
    n = 82

    scratch_thickness = .8
    
    img = "jupiter2.png"
    # load the image
    image = Image.open(img)
    # convert image to numpy array
    img_data = np.asarray(image)
    img_shape = img_data.shape
    img_scratch = np.array(img_data)
    scratch = np.zeros(img_shape, dtype=np.uint8)

    data = [int(i) for i in list(img_data.flatten())]

    # save dim of image
    y_max, x_max = image.size
    ex_max = math.ceil(3*x_max/k)*n
    ey_max = math.ceil(3*y_max/k)*n

    grs = ReedSolomon.GeneralizedReedSolomon(f=f, k=k, n=n, conventional_creation=True,
                                             alpha=0x2)

    # Encode image
    # We want to encode groups of 24 bytes, or 8 pixels and we will encode these to 30 bytes
    chunks = split(data, k)
    encoded_msg = []
    for i in tqdm(range(len(chunks))):
        encoded_msg += grs.encode(chunks[i])
    # Create a random scratch in the form ax + b
    a = random.randint(-1*3, 3)
    b = random.randint(0, y_max)
    a = -1
    b = 500

    # Something is weird with this fix that
    for i in range(0, len(encoded_msg), 3):
        x, y = divmod(i, ey_max)  # x - row, y - column
        y = int(k*(y/n)/3)
        if abs(y - (a*x+b)) != 0:
            prob = scratch_thickness*(y_max - abs(y - (a*x+b)))/(7*abs(y - (a*x+b))*abs(y - (a*x+b)) * y_max)
        else:
            if scratch_thickness != 0:
                prob = 1
            else:
                prob = 0
        if random.random() < prob and x < x_max and y < y_max:
            img_scratch[x, y] = np.array([255, 255, 255])
            scratch[x, y] = np.array([255, 255, 255])
            encoded_msg[i] = 0
            encoded_msg[i + 1] = 0
            encoded_msg[i + 2] = 0

    Image.fromarray(img_scratch).save('img_scratch.png')
    Image.fromarray(scratch).save('scratch.png')

    # Decode
    encoded_chunks = split(encoded_msg, n)
    decoded_msg = []
    for i in tqdm(range(len(encoded_chunks))):
        decoded_msg += grs.decode(encoded_chunks[i])

    decoded_img = split(split(split(decoded_msg, 3), img_shape[1]), img_shape[0])
    Image.fromarray(np.array(decoded_img[0], dtype=np.uint8)).save('img_decoded.png')


