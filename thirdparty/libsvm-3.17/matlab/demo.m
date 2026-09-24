gamma = 0.5;
K1 = Kernel(gamma,traindata(1:17,:),traindata(1:17,:));
n = size(K1,1);
K1 = [(1:n)',K1];
model = svmtrain(trainGID(1:17),K1,'-t 4');

K2 = Kernel(gamma,testdata(1:26,:),traindata(1:17,:));
n = size(K2,1);
K2 = [(1:n)',K2];
[predict_label, accuracy, dec_values] = svmpredict(testGID(1:26), K2, model);

[C,IA,~] = unique(testGID(1:26),'rows','R2012a');
k = size(C,1);
O = predict_label;
Y = testGID(1:26);
for i=1:k

   for j=1:k

      confusion_matrix(i,j) = length(find(Y == i & O == j));      

   end

end