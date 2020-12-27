import DNN_2layer
import DNN_3layer
import DNN_4layer


# test if the code works
if __name__ == '__main__':

    layer1 = layer2 = layer3 = layer4 = [10]
    DNN_2layer.DL_tune(layer1, layer2)
    DNN_3layer.DL_tune(layer1, layer2, layer3)
    DNN_4layer.DL_tune(layer1, layer2, layer3, layer4)