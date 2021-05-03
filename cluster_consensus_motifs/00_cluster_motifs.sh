#generate motif similarity file with HOCOMOCOv11_core_MOUSE.all.txt
sh run_MOUSE_tomtom.sh
#generate affinity propagation clustered motifs
python create-consensus-tomtom.py HOCOMOCOv11_core_MOUSE.all.txt HOCOMOCOv11_core_pcms_MOUSE_mono.txt HOCOMOCOv11_core_MOUSE_cluster
