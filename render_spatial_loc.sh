set -o xtrace
ORDER=(FID3DSplit FID3D FID3DGrads FID3DGradsJoinSig)
parallel=5
for i in 0 1 
do
    scenes=("${ORDER[@]:i:parallel}")
    echo $scenes
    for scene in ${scenes[@]}
    do
        manim -qh --renderer=opengl --write_to_movie spins.py $scene & 
    done
done
wait
manim-slides convert --to=pptx -cwidth=600 -cheight=600 ${ORDER[@]} FIDGrad.pptx
