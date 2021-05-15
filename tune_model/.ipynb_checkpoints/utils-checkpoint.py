import tensorflow as tf
import pandas as pd
import datetime
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score

def load_data(compdata_path):
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f') + ": Loading data ...")
    papers = pd.read_csv(compdata_path, index_col=False)
    papers.drop(['Unnamed: 0'], axis=1, inplace=True)

    return papers.fillna(0)


def create_model(n_features, layers):
    # create model
    model = tf.keras.Sequential()

    # add first hidden layer (and input layer)
    model.add(tf.keras.layers.Dense(layers[0], input_dim=n_features, activation='relu'))

    # add other hidden layers
    for i in range(1, len(layers)):
        model.add(tf.keras.layers.Dense(layers[i], activation='relu'))

    # add output layer
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))

    # Compile model
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=["accuracy"])
    return model


def dl_modeling(X_train, X_test, y_train, y_test, features, layers, epochs=10, batch_size=64):
    model = create_model(len(features), layers)
    md = model.fit(X_train[features], y_train, epochs=epochs, batch_size=64)

    pred = model.predict(X_test[features])

    fpr_dl, tpr_dl, _ = roc_curve(y_test, pred)
    return auc(fpr_dl, tpr_dl), accuracy_score(y_test, pred>0.5), precision_score(y_test, pred>0.5),\
           recall_score(y_test, pred>0.5)