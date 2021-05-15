# test run
folder_path='/Users/xuxiaopeng/Projects/data-covid-review/2021-01-25/'
test=1
epochs_list="1"
echo $epochs_list

# physical run
#folder_path='/home/xiaopeng/covid-archive/2021-02-07/'
#test=0
#epochs_list="5 10 15 20 30"
now=$(date +"%T")
echo "$now DL model parameter tuning started!"

for epochs in $epochs_list;
  do
    echo "Params: \n\tFolder path: $folder_path \n\tTest: $test, epochs: $epochs."
    python DNN_tune.py --epochs $epochs --hidden_layer_num 2 --test $test $folder_path 2>&1 >tune_2l_model_$epochs.out &
    python DNN_tune.py --epochs $epochs --hidden_layer_num 3 --test $test $folder_path 2>&1 >tune_3l_model_$epochs.out &
    python DNN_tune.py --epochs $epochs --hidden_layer_num 4 --test $test $folder_path 2>&1 >tune_4l_model_$epochs.out &
  done

wait
now=$(date +"%T")
echo "$now DL model parameter tuning finished!"
