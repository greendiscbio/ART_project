print("H0: Autoencoder accuracy results are higher than Logistic Regresion acuracy results.")
print("H1: Logistic Regresion accuracy results are higher than Autoencoder acuracy results.")

autoencoder_input =      [1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0]
autoencoder_prediction = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 
lr_input =     [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1]
lr_prediction = [1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1]

cont = 0
for i in range (len(autoencoder_input)):
    if(autoencoder_input[i]== autoencoder_prediction[i] and lr_input[i] != lr_prediction[i]):
        cont = cont+1

print('p-value H0 = ' + str(cont/len(autoencoder_input)))

cont = 0
for i in range (len(autoencoder_input)):
    if(autoencoder_input[i]!= autoencoder_prediction[i] and lr_input[i] == lr_prediction[i]):
        cont = cont+1

print('p-value H1 = ' + str(cont/len(autoencoder_input)))