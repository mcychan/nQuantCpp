# nQuantCpp
nQuantCpp includes top 4 color quantizers producing high quality 256 color 8 bit PNG images using an algorithm optimized for the highest quality possible.

Pairwise Nearest Neighbor quantization, 
NeuQuant Neural-Net Quantization Algorithm, 
Xialoin Wu's fast optimal color quantizer, 
DL3 Quantization

Each quantization algorithm has its own advantages. Pairwise Nearest Neighbor quantization minimized color loss for photo having red lips and supports 256 or less colors. NeuQuant Neural-Net Quantization Algorithm produces smooth photo quantization especially for natual landscape photo. Only Xialoin Wu's fast optimal color quantizer fully support image having transparent color. DL3 Quantization supports 256 or less colors. nQuantCpp also provides a command line wrapper in case you want to use it from the command line.

Either download nQuant from this site or add it to your Visual Studio project seamlessly.
The main features show up to discuss would be the error diffusion and dithering.

Here some examples of output:

Original image<img src="https://i.stack.imgur.com/fOcIL.png" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/4jLkg.png" /><br>
Reduced to 256 colors by Xialoin Wu's fast optimal color Quantization Algorithm<img src="https://i.stack.imgur.com/PGWVF.png" /><br><br>
Original photo<img src="https://i.stack.imgur.com/jFvEG.jpg" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/yJIjQ.gif" />
Reduced to 256 colors by Pairwise Nearest Neighbor Quantization Algorithm<img src="https://i.stack.imgur.com/dPTml.gif" />

Welcome for C++ experts for further improvement or provide color quantization algorithms better than the above algorithms.
