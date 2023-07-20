module load METAL/11.3.25
for d in */ ; do
	cd "$d";
	metal source ../metal_script.txt > metal_${d%/}.log;
	cd ..;
done
