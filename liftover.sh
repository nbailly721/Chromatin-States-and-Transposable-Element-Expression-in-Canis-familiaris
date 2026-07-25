!/bin/bash
module load kent_tools/486

liftOver OV_13_dense.bed canFam3ToGCF_014441545.1.over.chain OV_ROS.bed OV_unmapped.bed

liftOver SP_13_dense.bed canFam3ToGCF_014441545.1.over.chain SP_ROS.bed SP_unmapped.bed

liftOver CR_13_dense.bed canFam3ToGCF_014441545.1.over.chain CR_ROS.bed CR_unmapped.bed

liftOver CL_13_dense.bed canFam3ToGCF_014441545.1.over.chain CL_ROS.bed CL_unmapped.bed
