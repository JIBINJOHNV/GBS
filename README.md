# GBS


python GBS_Double_Digestion.py \
-FastaPath /data/NGC_Data/Analysis/NGC_Internal/GBSDATA_21022022/Mouse_inslico/IDX_Fasta/PipelineAutomation/ \
-FastFile genome.fas \
-R1Site G^AATTC -R2Site T^TAA \
-R1Enzyme EcoRI -R2Enzyme MspI \
-LSize 90 -USize 130 \
-MQuality 20 -OutName Mouse
