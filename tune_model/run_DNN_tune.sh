# Deep learning model tuning
folder_path='/home/xiaopengxu/Desktop/data-covid-review/2021-05-11/'
test=0 # 1 for test run, 0 for actually tuning
epochs_list="5 10"
echo $epochs_list

now=$(date +"%T")
echo "$now DL model parameter tuning started!"

for epochs in $epochs_list;
  do
    echo "Params:  folder path: $folder_path, test: $test, epochs: $epochs."
    python DNN_tune.py --epochs $epochs --hidden_layer_num 2 --test $test $folder_path 2>&1 >tune_2l_model_$epochs.out &
    python DNN_tune.py --epochs $epochs --hidden_layer_num 3 --test $test $folder_path 2>&1 >tune_3l_model_$epochs.out &
    python DNN_tune.py --epochs $epochs --hidden_layer_num 4 --test $test $folder_path 2>&1 >tune_4l_model_$epochs.out &
  done

wait

now=$(date +"%T")
echo "$now DL model parameter tuning finished!"
