import os
from PIL import Image, ImageDraw, ImageColor

line_range = (402, 745)
colour_range = (320, -40)

SATURATION = 100
VALUE = 100

lines = [419.8, 426.7, 433.3, 449.8, 457.2, 468.6, 518.4, 544.9, 559.6, 592.0, 602.0, 695.1]

HEIGHT = 100
SCALE = 4


with Image.new("HSV", ((line_range[1] - line_range[0])*SCALE, HEIGHT)) as img:
    canvas = ImageDraw.Draw(img)

    for line in lines:
        hue = int((((line-line_range[0]) / (line_range[1]-line_range[0])) * (colour_range[1]-colour_range[0]) + colour_range[0]) * (255/360))
        canvas.line([((line-line_range[0])*SCALE, 0), ((line-line_range[0])*SCALE, HEIGHT)], width=SCALE, fill=(hue, 255, 255))
    
    img.show()
    img = img.convert("RGB")

    img.save("spectrum.jpg", "JPEG")