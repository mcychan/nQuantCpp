# nQuantCpp
nQuantCpp includes top 6 color quantization algorithms for visual c++ producing high quality optimized images. I enhance each of the algorithms to support semi transparent images. 
nQuantCpp also provides a command line wrapper in case you want to use it from the command line.

Either download nQuantCpp from this site or add it to your Visual Studio project seamlessly.
PNG is useful because it's the only widely supported format which can store partially transparent images. The format uses compression, but the files can still be large. Use Color quantization algorithms can be chosen by command line since version 1.10 using the /a algorithm.
Only png can support semi transparent image and desired color depth. Gif can ensure the number of colors for the converted image is 256 or less. Bmp does support desired color depth. Jpg only supports 24-bit image format.

Let's climb up the mountain: Ready, Go!!!

<p>Original photo of climbing<br /><img src="https://mcychan.github.io/PnnQuant.js/demo/img/climb.jpg" /></p>
<p>Reduced to 256 colors by Divisive hierarchical clustering algorithm<br /><img src="https://i.stack.imgur.com/Qitc4.png" /></p>
<p>Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<br /><img src="https://i.stack.imgur.com/ebUOv.png" /></p>
<p>Reduced to 16 colors by Fast pairwise nearest neighbor based algorithm<br /><img src="https://i.stack.imgur.com/07EFv.png" /></p>
<p>Reduced to 16 colors by Fast pairwise nearest neighbor based algorithm with CIELAB color space<br /><img src="https://i.stack.imgur.com/6GxLY.png" /></p>
<p>Reduced to 16 colors by Xialoin Wu's fast optimal color Quantization Algorithm<br /><img src="https://i.stack.imgur.com/8PEDu.png" /></p>
<hr />
<p>Original photo of Aetna's Hartford headquarters<br /><img src="https://mcychan.github.io/PnnQuant.js/demo/img/SE5x9.jpg" /></p>
<p>Reduced to 256 colors by NeuQuant Neural-Net Quantization Algorithm<br /><img src="https://i.stack.imgur.com/0sDDn.png" /></p>
<p>Reduced to 256 colors by Fast pairwise nearest neighbor based algorithm<br /><img src="https://i.stack.imgur.com/SB6NJ.png" /></p><hr>

<p>Original image of Hong Kong Cuisines<br /><img src="https://mcychan.github.io/PnnQuant.js/demo/img/old-HK.jpg" /></p>
<b><a href="http://www.cs.joensuu.fi/sipu/pub/Threshold-JEI.pdf">Fast pairwise nearest neighbor based algorithm</a></b> with CIELAB color space in 16 colors<br>
High quality and fast<br />
<img src="https://repository-images.githubusercontent.com/121180544/95f03a00-7780-11eb-9330-f2e6232c283f" alt="Fast pairwise nearest neighbor based algorithm with CIELAB color space with CIELAB color space in 16 colors"></p>
<p><b><a href="http://cg.cs.tsinghua.edu.cn/people/~huanghz/publications/TIP-2015-CombinedColorQuantization.pdf">Efficient, Edge-Aware, Combined Color Quantization and Dithering </a></b> with CIELAB color space in 16 colors<br />
Higher quality for 32 or less colors but slower<br />
<img src="https://i.stack.imgur.com/Vk0dG.png" alt="Efficient, Edge-Aware, Combined Color Quantization and Dithering with CIELAB color space in 16 colors"></p>
<p><b><a href="https://people.eecs.berkeley.edu/~dcoetzee/downloads/scolorq/">Spatial color quantization</a></b> with CIELAB color space in 16 colors<br />
Higher quality for 32 or less colors but the slowest<br />
<img src="https://i.stack.imgur.com/YB3hZ.png" alt="Spatial color quantization with CIELAB color space in 16 colors"></p>
All in all, the top 3 color quantization algorithms for 256 colors are:
<ol>
<li>Fast pairwise nearest neighbor based algorithm</li>
<li>Xialoin Wu's fast optimal color quantizer</li>
<li>NeuQuant Neural-Net Quantization Algorithm with CIELAB color space</li>
</ol>
The top 3 color quantization algorithms for 32 colors or less are:
<ol>
<li>Fast pairwise nearest neighbor based algorithm with CIELAB color space advanced version</li>
<li>Efficient, Edge-Aware, Combined Color Quantization and Dithering algorithm with CIELAB color space</li>
<li>Spatial color quantization algorithm with CIELAB color space</li>
</ol>

If you are using the command line. Assuming you are in the same directory as nQuantCpp.exe, you would enter: nQuantCpp yourImage.jpg /m 16
To avoid dot gain, /d n can set the dithering to false. However, false contours will be resulted for gradient color zones.
nQuantCpp will quantize yourImage.jpg and create yourImage-PNNLABquant16.png in the same directory.

The readers can see coding of the error diffusion and dithering are quite similar among the above quantization algorithms. 
Each algorithm has its own advantages. I share the source of color quantization to invite further discussion and improvements.
Such source code are written in C++ to gain best performance. It is readable and convertible to <a href="https://github.com/mcychan/nQuant.cs">c#</a>, <a href="https://github.com/mcychan/nQuant.j2se">java</a>, or <a href="https://github.com/mcychan/PnnQuant.js">javascript</a>.
Welcome for C++ experts for further improvement or provide color quantization algorithms better than the above algorithms.
Please use issues to track ideas, enhancements, tasks, or bugs for work on GitHub.
