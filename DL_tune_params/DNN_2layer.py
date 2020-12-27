import numpy as np
import tensorflow as tf
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc


def create_model(n_features, layers):
    # create model
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Dense(layers[0], input_dim=n_features, activation='relu'))

    # add hidden layers
    for l in range(1, len(layers)):
        model.add(tf.keras.layers.Dense(50, activation='relu'))
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))

    # Compile model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=["accuracy"])
    return model


def DL_Modeling(features, layers, epochs):
    X = pd.read_csv('../data/lda_doc2vec_feats.csv', index_col=False)
    y = pd.read_csv('../data/lda_doc2vec_targs.csv', index_col=False)
    X = X.drop(columns=['Unnamed: 0'])
    y = y.drop(columns=['Unnamed: 0'])
    X = X.fillna(0)

    df_emb = pd.read_csv('../data/graph_embedding.csv', index_col=False)
    df_emb = df_emb.drop(columns=['Unnamed: 0'])
    X = pd.concat([X, df_emb], axis=1)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    model = create_model(len(features), layers)
    md = model.fit(X_train[features], y_train, epochs=epochs, batch_size=64)
    pred = model.predict(X_test[features])

    fpr_dl, tpr_dl, _ = roc_curve(y_test, pred)

    return fpr_dl, tpr_dl


def DL_tune(layer1, layer2):
    em_size = 50
    F_Doc2Vec = [str(i) for i in range(em_size)]
    F_gem = ['em ' + str(i + 1) for i in range(em_size)]

    print('2 layer tuning starts: total tests {0}'.format(str(len(layer1) * len(layer2))))
    file1 = open("../data/2layer_aucs.csv", "w")
    file1.write("layer1\tlayer2\tauc\n")
    for i2 in range(len(layer2)):
        for i1 in range(len(layer1)):
            fpr_dl, tpr_dl = DL_Modeling(F_Doc2Vec + F_gem, [layer1[i1], layer2[i2]], 1)
            file1.write("{0}\t{1}\t{2}\n".format(str(layer1[i1]), str(layer2[i2]), str(
                auc(fpr_dl, tpr_dl))))

    file1.close()
    print('2 layer tuning finished!')


# Two layer params tuning
if __name__ == '__main__':
    layer1 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140]
    layer2 = [10, 20, 30, 40, 50, 60, 70, 80]

    DL_tune(layer1, layer2)

