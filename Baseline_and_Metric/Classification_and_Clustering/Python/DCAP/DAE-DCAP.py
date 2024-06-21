

import tensorflow as tf
import numpy as np  
import matplotlib.pyplot as plt  
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

with open(r"simulation.csv", 'r') as f:
    data = pd.read_csv(f)

#print(data.shape)
tcga_input=np.transpose(data)
print(tcga_input.shape[1])
length1 = tcga_input.shape[1]
learning_rate = 0.0001
training_epochs = 100
batch_size = 125
display_step = 2
examples_to_show = 10
n_input = tcga_input.shape[1]
  

X = tf.placeholder("float", [None, n_input])  
  

n_hidden_1 = 200
n_hidden_2 = 50
n_hidden_3 = 2

weights = {  
    'encoder_h1': tf.Variable(tf.random_normal([n_input, n_hidden_1])),  
    'encoder_h2': tf.Variable(tf.random_normal([n_hidden_1, n_hidden_2])),
    'encoder_h3': tf.Variable(tf.random_normal([n_hidden_2, n_hidden_3])),
    'decoder_h1': tf.Variable(tf.random_normal([n_hidden_3, n_hidden_2])),
    'decoder_h2': tf.Variable(tf.random_normal([n_hidden_2, n_hidden_1])),
    'decoder_h3': tf.Variable(tf.random_normal([n_hidden_1, n_input])),


}  
biases = {  
    'encoder_b1': tf.Variable(tf.random_normal([n_hidden_1])),  
    'encoder_b2': tf.Variable(tf.random_normal([n_hidden_2])),
    'encoder_b3': tf.Variable(tf.random_normal([n_hidden_3])),
    'decoder_b1': tf.Variable(tf.random_normal([n_hidden_2])),
    'decoder_b2': tf.Variable(tf.random_normal([n_hidden_1])),
    'decoder_b3': tf.Variable(tf.random_normal([n_input])),
}  
  

def encoder(x):  
    layer_1 = tf.nn.tanh(tf.add(tf.matmul(x, weights['encoder_h1']),
                                   biases['encoder_b1']))  
    layer_2 = tf.nn.tanh(tf.add(tf.matmul(layer_1, weights['encoder_h2']),
                                   biases['encoder_b2']))
    layer_3 = tf.nn.tanh(tf.add(tf.matmul(layer_2, weights['encoder_h3']),
                                   biases['encoder_b3']))
    return layer_3
  
  

def decoder(x):  
    layer_1 = tf.nn.tanh(tf.add(tf.matmul(x, weights['decoder_h1']),
                                   biases['decoder_b1']))  
    layer_2 = tf.nn.tanh(tf.add(tf.matmul(layer_1, weights['decoder_h2']),
                                   biases['decoder_b2']))
    layer_3 = tf.nn.tanh(tf.add(tf.matmul(layer_2, weights['decoder_h3']),
                                   biases['decoder_b3']))
    return layer_3
  

encoder_op = encoder(X)  
decoder_op = decoder(encoder_op)  
  

y_pred = decoder_op  
y_true = X  
  

cost = tf.reduce_mean(tf.pow(y_true - y_pred, 2))
optimizer = tf.train.AdamOptimizer(learning_rate).minimize(cost)  
  
with tf.Session() as sess:
    if int((tf.__version__).split('.')[1]) < 12 and int((tf.__version__).split('.')[0]) < 1:  
        init = tf.initialize_all_variables()  
    else:  
        init = tf.global_variables_initializer()  
    sess.run(init)  

    total_batch = int(len(tcga_input)/batch_size)
    for epoch in range(training_epochs):  
        for i in range(total_batch):
            batch_xs = tcga_input[((i) * batch_size):((i + 1) * batch_size)] + 0.3 * np.random.rand(length1)   #added nosie
            # Run optimization op (backprop) and cost op (to get loss value)
            _, c = sess.run([optimizer, cost], feed_dict={X: batch_xs})
        if epoch % display_step == 0:  
            print("Epoch:", '%04d' % (epoch+1), "cost=", "{:.9f}".format(c))
        if epoch == training_epochs - 1:
                fea_output = sess.run([encoder_op], feed_dict={X: tcga_input})
                # print(fea_output)
                print(np.array(fea_output).shape)
                np.savetxt(r'fea.csv', np.array(fea_output[0]), delimiter=',')
                dd = np.array(fea_output[0])
    print("Optimization Finished!")
    print(dd.shape)
    clf = KMeans(n_clusters=2)
    clf.fit(dd)
    centers = clf.cluster_centers_
    labels = clf.labels_
    silhouetteScore = silhouette_score(dd, labels, metric='euclidean')
    print(centers)
    print(silhouetteScore)
    # encode_decode = sess.run(
    #     y_pred, feed_dict={X: mnist.test.images[:examples_to_show]})
    # f, a = plt.subplots(2, 10, figsize=(10, 2))  #return fig，axes
    # for i in range(examples_to_show):  
    #     a[0][i].imshow(np.reshape(mnist.test.images[i], (28, 28)))  
    #     a[1][i].imshow(np.reshape(encode_decode[i], (28, 28)))  
    # plt.show() 