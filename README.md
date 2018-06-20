# nQuantCpp
nQuantCpp includes top 5 color quantizers producing high quality 256 color 8 bit PNG images using an algorithm optimized for the highest quality possible.

Fast pairwise nearest neighbor based algorithm, 
NeuQuant Neural-Net Quantization Algorithm, 
Xialoin Wu's fast optimal color quantizer, 
DL3 Quantization,
Spatial color quantization

Fast pairwise nearest neighbor based algorithm minimized color loss for photo having red lips and supports 256 or less colors with transparency. NeuQuant Neural-Net Quantization Algorithm produces smooth photo quantization especially for natual landscape photo and fully supports image having transparent color. Xialoin Wu's fast optimal color quantizer fully supports image having transparent color. DL3 Quantization supports 256 or less colors. Spatial color quantization supports 64 or less colors. nQuantCpp also provides a command line wrapper in case you want to use it from the command line.

Either download nQuant from this site or add it to your Visual Studio project seamlessly.
The main features show up to discuss would be the error diffusion and dithering.

Here some examples of output:

Original image<img src="https://i.stack.imgur.com/fOcIL.png" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/4jLkg.png" /><br>
Reduced to 256 colors by Xialoin Wu's fast optimal color Quantization Algorithm<img src="https://i.stack.imgur.com/PGWVF.png" /><br><br>
Original photo<img src="https://i.stack.imgur.com/jFvEG.jpg" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/yJIjQ.gif" />
Reduced to 256 colors by Fast pairwise nearest neighbor based algorithm<img src="https://i.stack.imgur.com/dPTml.gif" />

The readers can see coding of the error diffusion and dithering are quite similar among the above quantization algorithms. 
Each algorithm has its own advantages. I share the source of color quantization to invite further discussion and improvements.
Such source code are written in C++ to gain best performance. It is readable and convertible to c#, java, or javascript.
Welcome for C++ experts for further improvement or provide color quantization algorithms better than the above algorithms.
Please send email to miller.chan@gmail.com to give suggestions or request the advanced version of Pairwise Nearest Neighbor quantization.
