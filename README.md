# nQuantCpp
nQuant cpp includes top 4 color quantizers producing high quality 256 color 8 bit PNG images using an algorithm optimized for the highest quality possible.

Pairwise Nearest Neighbor quantization, 
NeuQuant Neural-Net Quantization Algorithm, 
Xialoin Wu's fast optimal color quantizer, 
DL3 Quantization

Using nQuant cpp
Each quantization pros and cons. Pairwise Nearest Neighbor quantization no color loss for photo having red lips and support 256 or less colors. NeuQuant Neural-Net Quantization Algorithm offer smooth photo quantization especially for natual landscape photo. Only Xialoin Wu's fast optimal color quantizer fully support transparency. DL3 Quantization supports 256 or less colors. nQuantCpp also provides a command line wrapper in case you want to use it from the command line. To get started:

Either download nQuant from this site or add it to your Visual Studio project seamlessly.
Then add nQuant to your project and add a reference to it.
If you are using native C++, you would call nQuant cpp as follows:
 

 auto m_pImage256Color = make_unique<Bitmap>(w, h, PixelFormat8bppIndexed);
 auto m_pImage = ConvertTo(m_pImage.get(), PixelFormat32bppARGB);

 bool bSucceeded = wuQuantizer.QuantizeImage(m_pImage.get(), m_pImage256Color.get(), 256);

Welcome for C++ experts for further improvement or provide color quantization algorithms better than the above algorithms.
