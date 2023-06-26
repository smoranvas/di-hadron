#int nEvent = atoi(argv[1]);
#int Z = atoi(argv[2]);
#int A = atoi(argv[3]);
#int mode = atof(argv[4]);
#double K = atof(argv[5]);
#std::string path(argv[6]);
#auto header = std::string(argv[7]);
#auto configfile = std::string(argv[8]);

nEvent=1000

K=4
table_path=table/
mode=1
config_file=eA_5GeV.txt


#$exe=../../eHIJING-pythia/eHIJING-examples/build/clas
exe=build/clas_dihadron

Z=1
A=2
header=test_D
mkdir $header
$exe $nEvent $Z $A $mode $K $table_path $header $config_file


Z=6
A=12
header=test_C
mkdir $header
$exe $nEvent $Z $A $mode $K $table_path $header $config_file

Z=26
A=56
header=test_Fe
$exe $nEvent $Z $A $mode $K $table_path $header $config_file

Z=82
A=207
header=test_Fe
$exe $nEvent $Z $A $mode $K $table_path $header $config_file
