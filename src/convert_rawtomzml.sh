aa=$(pwd)

pwd=$(pwd)
#export ex_dir='/Users/hock-sengnguan/Documents/research/sugar_new/programs/others/test1/ThermoRawFileParser1.4.2'
export ex_dir='./ThermoRawFileParser1.4.2'
export raw_dir=$pwd
export out_dir=$pwd

mono $ex_dir/ThermoRawFileParser.exe -i=$raw_dir/$1 -o=$out_dir/ -f=1 -p=1- -m=1

#mono $ex_dir/ThermoRawFileParser.exe query -i=$raw_dir/$1 -b=$out_dir/$2 -p  -n=$3 -l=0

