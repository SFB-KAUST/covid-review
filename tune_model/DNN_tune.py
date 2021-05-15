import argparse
import datetime
from sklearn.model_selection import train_test_split
from utils import load_data, dl_modeling


def tune_2layer(X_train, X_test, y_train, y_test,feature_set, save_path, layer1, layer2, epochs=10):
    """
        Two layer params tuning
    :param papers: pandas dataframe. Input training data
    :param save_path: str. Path to save tuning results
    :param layer1: list. Number of nodes in the first layer
    :param layer2: list. Number of nodes in the second layer
    :param epochs: int. Tensorflow epochs, default 10
    :return: No return value
    """

    for feature in feature_set.keys():
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
              + ":Start tuning 2 hidden layer model with {}".format(feature)
              + ": total runs {0}".format(len(layer1) * len(layer2)))

        file1 = open(save_path + "tune_dl_parmas_{}.2_hid_layers_{}.csv".format(feature, epochs), "w")
        file1.write("layer1\tlayer2\tauc\tacc\tprec\trecall\n")
        for i2 in range(len(layer2)):
            for i1 in range(len(layer1)):
                auc, acc, prec, recall = dl_modeling(X_train, X_test, y_train, y_test, feature_set[feature],
                                                     [layer1[i1], layer2[i2]], epochs=epochs)
                file1.write("{0}\t{1}\t".format(layer1[i1], layer2[i2])
                            + "{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(auc, acc, prec, recall))
        file1.close()
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + ' 2 layer tuning finished!')


def tune_3layer(X_train, X_test, y_train, y_test, feature_set, save_path, layer1, layer2, layer3, epochs=10):
    """
        3 layer params tuning
    :param papers: pandas dataframe. Input training data
    :param save_path: str. Path to save tuning results
    :param layer1: list. Number of nodes in the 1st layer
    :param layer2: list. Number of nodes in the 2nd layer
    :param layer3: list. Number of nodes in the 3rd layer
    :param epochs: int. Tensorflow epochs, default 10
    :return: No return value
    """

    for feature in feature_set.keys():
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
              + ":Start tuning 3 hidden layer model with {}".format(feature)
              + ": total runs {0}".format(len(layer1) * len(layer2) * len(layer3)))

        file1 = open(save_path + "tune_dl_parmas_{}.3_hid_layers_{}.csv".format(feature, epochs), "w")
        file1.write("layer1\tlayer2\tlayer3\tauc\tacc\tprec\trecall\n")
        for i3 in range(len(layer3)):
            for i2 in range(len(layer2)):
                for i1 in range(len(layer1)):
                    auc, acc, prec, recall = dl_modeling(X_train, X_test, y_train, y_test, feature_set[feature],
                                                         [layer1[i1], layer2[i2], layer3[i3]], epochs=epochs)
                    file1.write("{0}\t{1}\t{2}\t".format(layer1[i1], layer2[i2], layer3[i3])
                                + "{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(auc, acc, prec, recall))
        file1.close()
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + ' 3 layer tuning finished!')


def tune_4layer(X_train, X_test, y_train, y_test, feature_set, save_path, layer1, layer2, layer3, layer4, epochs=10):
    """
        4 layer params tuning
    :param papers: pandas dataframe. Input training data
    :param save_path: str. Path to save tuning results
    :param layer1: list. Number of nodes in the 1st layer
    :param layer2: list. Number of nodes in the 2nd layer
    :param layer3: list. Number of nodes in the 3rd layer
    :param layer4: list. Number of nodes in the 4th layer
    :param epochs: int. Tensorflow epochs, default 10
    :return: No return value
    """

    for feature in feature_set.keys():
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
              + ":Start tuning 4 hidden layer model with {}".format(feature)
              + ": total runs {0}".format(len(layer1) * len(layer2) * len(layer3) * len(layer4)))

        file1 = open(save_path + "tune_dl_parmas_{}.4_hid_layers_{}.csv".format(feature, epochs), "w")
        file1.write("layer1\tlayer2\tlayer3\tlayer4\tauc\tacc\tprec\trecall\n")
        for i4 in range(len(layer4)):
            for i3 in range(len(layer3)):
                for i2 in range(len(layer2)):
                    for i1 in range(len(layer1)):
                        auc, acc, prec, recall = dl_modeling(X_train, X_test, y_train, y_test, feature_set[feature],
                                                             [layer1[i1], layer2[i2], layer3[i3], layer4[i4]],
                                                             epochs=epochs)
                        file1.write("{0}\t{1}\t{2}\t{3}\t".format(layer1[i1], layer2[i2], layer3[i3], layer4[i4])
                                    + "{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(auc, acc, prec, recall))
        file1.close()
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + ' 4 layer tuning finished!')


def main(args):
    layer1 = layer2 = layer3 = layer4 = [10]
    if args.test != 1:                                                           # if not test
        layer1 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
        layer2 = [10, 20, 30, 40, 50, 60, 70, 80]
        layer3 = [10, 20, 30, 40, 50]
        layer4 = [10, 20, 30]

    inoutpath = args.folder_path
    compdata_path = inoutpath + 'features.ori_doc2vec_deepwalk_scores_train.csv'

    X = load_data(compdata_path)
    y = X.is_published

    em_size = 50
    F_Doc2Vec = ['dv ' + str(i + 1) for i in range(em_size)]
    F_DeepWalk = ['dw ' + str(i + 1) for i in range(em_size)]
    F_ascore = ['auth.papers.score', 'auth.citations.score', 'art.citations.score', 'art.views.score', 'final.score']
    
    feature_set = {"Doc2Vec": F_Doc2Vec,
                   "DeepWalk": F_DeepWalk,
                   "Ascore": F_ascore,
                   "Doc2Vec+DeepWalk+AScore": F_Doc2Vec + F_DeepWalk + F_ascore}
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    
    if args.hidden_layer_num == 2:
        tune_2layer(X_train, X_test, y_train, y_test, feature_set, inoutpath, layer1, layer2, epochs=args.epochs)
    elif args.hidden_layer_num == 3:
        tune_3layer(X_train, X_test, y_train, y_test, feature_set, inoutpath, layer1, layer2, layer3, epochs=args.epochs)
    elif args.hidden_layer_num == 4:
        tune_4layer(X_train, X_test, y_train, y_test, feature_set, inoutpath, layer1, layer2, layer3, layer4, epochs=args.epochs)

    return


def parse_args():
    parser = argparse.ArgumentParser(description="Python script to tune DL model layers")
    parser.add_argument("folder_path", type=str, help="Path to the folder containing training data")
    parser.add_argument("--test", type=int, default=1, help="Test mode if 1, default 1")
    parser.add_argument("--hidden_layer_num", type=int, default=2,
                        help="Number of hidden layers, in {2,3,4}, default 2")
    parser.add_argument("--epochs", type=int, default=10, help="Epochs in Tensorflow, default 10")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)



