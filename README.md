# nQuantCpp
nQuantCpp includes top 6 color quantization algorithms for visual c++ producing high quality optimized images. I enhance each of the algorithms to support semi transparency images. 

Divisive hierarchical clustering algorithm utilizes the commonly used binary splitting strategy along with several carefully selected heuristics that ensure a good balance between effectiveness and efficiency.
Fast pairwise nearest neighbor based algorithm is log linear and considerably faster than optimal thresholding, it is applicable in real-time image processing applications. 
NeuQuant Neural-Net Quantization Algorithm is a self-organizing Kohonen neural network for quantizing colour graphics images.
Xialoin Wu's fast optimal color quantizer, which is one of the most effective color quantization methods, provides excellent results.
Efficient, Edge-Aware, Combined Color Quantization and Dithering is a novel algorithm to simultaneously accomplish color quantization and dithering of images.
Spatial color quantization is a novel technique for palette selection and dithering with a simple perceptual model of human vision to produce superior results for many types of images.

Fast pairwise nearest neighbor based algorithm minimized color loss for photo having red lips and supports 256 or less colors. NeuQuant Neural-Net Quantization Algorithm produces smooth photo quantization especially for natual landscape photo supports image reduction to 64 or more colors. Xialoin Wu's fast optimal color quantizer supports image reduction to 64 or more colors. Efficient, Edge-Aware, Combined Color Quantization and Dithering, Spatial color quantization supports 64 or less colors. nQuantCpp also provides a command line wrapper in case you want to use it from the command line.

Either download nQuantCpp from this site or add it to your Visual Studio project seamlessly.
PNG is useful because it's the only widely supported format which can store partially transparent images. The format uses compression, but the files can still be large. Use Color quantization algorithms can be chosen by command line since version 1.10 using the /a algorithm.
Only png can support semi transparent image and desired color depth. Gif can ensure the number of colors for the converted image is 256 or less. Bmp does support desired color depth. Jpg only supports 24-bit image format.

Here some examples of output:

Original image<img src="http://i.imgur.com/h9ghTMB.png" /><br>
Reduced to 256 colors by Divisive hierarchical clustering algorithm<img src="https://i.stack.imgur.com/viRTI.png" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/G1jkp.png" /><br>
Reduced to 16 colors by Fast pairwise nearest neighbor based algorithm<img src="https://i.stack.imgur.com/ry1oi.png" /><br>
Reduced to 16 colors by Xialoin Wu's fast optimal color Quantization Algorithm<img src="https://i.stack.imgur.com/De9xw.png" /><br><br>
Original photo<br>
<img src="https://i.stack.imgur.com/SE5x9.png" /><br>
Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<img src="https://i.stack.imgur.com/0sDDn.png" /><br>
Reduced to 256 colors by Fast pairwise nearest neighbor based algorithm<img src="https://i.stack.imgur.com/SB6NJ.png" /><br>

<p>Original image<br><img src="https://i.stack.imgur.com/F90bn.jpg" /><br>
<b><a href="http://www.cs.joensuu.fi/sipu/pub/Threshold-JEI.pdf">Fast pairwise nearest neighbor based algorithm with CIELAB color space</a></b> with 16 colors<br>
High quality and fast<br>
<img src="https://i.stack.imgur.com/2kFxV.png" alt="Fast pairwise nearest neighbor based algorithm with CIELAB color space with 16 colors"></p>
<p><b><a href="http://cg.cs.tsinghua.edu.cn/people/~huanghz/publications/TIP-2015-CombinedColorQuantization.pdf">Efficient, Edge-Aware, Combined Color Quantization and Dithering</a></b> with 16 colors<br>
Higher quality for 32 or less colors but slower<br>
<img src="https://i.stack.imgur.com/cVYMP.png" alt="Efficient, Edge-Aware, Combined Color Quantization and Dithering with 16 colors"></p>
<p><b><a href="https://people.eecs.berkeley.edu/~dcoetzee/downloads/scolorq/">Spatial color quantization</a></b> with 16 colors<br>
Higher quality for 32 or less colors but the slowest<br>
<img src="https://i.stack.imgur.com/DVdGv.png" alt="Spatial color quantization with 16 colors"></p>

If you are using the command line. Assuming you are in the same directory as nQuantCpp.exe, you would enter: nQuantCpp yourImage.jpg /m 16

nQuantCpp will quantize yourImage.jpg and create yourImage-PNNLABquant16.png in the same directory.

The readers can see coding of the error diffusion and dithering are quite similar among the above quantization algorithms. 
Each algorithm has its own advantages. I share the source of color quantization to invite further discussion and improvements.
Such source code are written in C++ to gain best performance. It is readable and convertible to <a href="https://github.com/mcychan/nQuant.cs">c#</a>, <a href="https://github.com/mcychan/nQuant.j2se">java</a>, or <a href="https://github.com/mcychan/PnnQuant.js">javascript</a>.
Welcome for C++ experts for further improvement or provide color quantization algorithms better than the above algorithms.
Please send email to miller.chan@gmail.com to report issues or give suggestions.
