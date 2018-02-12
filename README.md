# nQuantCpp
nQuant cpp Color Quantizer
nQuant cpp is a visual c++ color quantizer producing high quality 256 color 8 bit PNG images using an algorithm optimized for the highest quality possible.

nQuant cpp was originally developed as part of nQuant .net which is an http module that automatically minifies and merges CSS as well as sprites their background images on the fly. I wanted the sprited files to be optimized and I was not satisfied with the size of the 32 bit images that the gdiplus library was producing nor was the quantization output of such quantizers as OctTree and WebSafe of acceptable quality. I re-encode this quantizer and the results are images 3x smaller than their 32 bit originals with practically no perceptible quality loss.

See this blog post for a description of this effort.

Why is nQuant's quality often better than other color quantizers?
nQuant is an adaptation of Xialoin Wu's fast optimal color quantizer see Graphics Gems vol. II, pp. 126-133. This algorithm often produces optimized images with less quality degradation than other quantizers due to its methodology of optimizing based on clusters of colors that are closely related to one another rather than simply finding the most used colors in an image. nQuant's algorithm takes Wu's three dimensional RGB based quantization strategy and adapts it to work with transparent PNGs.

Using nQuant cpp
Another advantage of nQuant is that it is a .net library that you can integrate nicely with your own .net code while many of the popular quantizers only provide command line implementations. nQuant also provides a command line wrapper in case you want to use it from the command line. To get started:

Either download nQuant from this site or add it to your Visual Studio project seamlessly.
Then add nQuant to your project and add a reference to it.
If you are using native C++, you would call nQuant cpp as follows:
 

 auto m_pImage256Color = make_unique<Bitmap>(w, h, PixelFormat8bppIndexed);
 auto m_pImage = ConvertTo(m_pImage.get(), PixelFormat32bppARGB);

 bool bSucceeded = wuQuantizer.QuantizeImage(m_pImage.get(), m_pImage256Color.get(), 256);
