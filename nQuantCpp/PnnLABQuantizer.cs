using System;
using System.Collections.Generic;
using System.Drawing;

namespace PnnQuant
{
    public class PnnLABQuantizer : PnnQuantizer
    {
        private readonly double PR = .299, PG = .587, PB = .114;
        private double ratio = 1.0;
        private readonly Dictionary<int, CIELABConvertor.Lab> pixelMap = new Dictionary<int, CIELABConvertor.Lab>();
        private sealed class Pnnbin
        {
            internal float ac, Lc, Ac, Bc;
            internal int cnt;
            internal int nn, fw, bk, tm, mtm;
            internal float err;
        }
        private void GetLab(int argb, out CIELABConvertor.Lab lab1)
        {
            if (!pixelMap.TryGetValue(argb, out lab1))
            {
                lab1 = CIELABConvertor.RGB2LAB(Color.FromArgb(argb));
                pixelMap[argb] = lab1;
            }
        }
        private void Find_nn(Pnnbin[] bins, int idx)
        {
            int nn = 0;
            var err = 1e100;

            var bin1 = bins[idx];
            int n1 = bin1.cnt;
            var lab1 = new CIELABConvertor.Lab
            {
                alpha = bin1.ac, L = bin1.Lc, A = bin1.Ac, B = bin1.Bc
            };
            for (int i = bin1.fw; i != 0; i = bins[i].fw)
            {
                float n2 = bins[i].cnt;
                float nerr2 = (n1 * n2) / (n1 + n2);
                if (nerr2 >= err)
                    continue;

                var lab2 = new CIELABConvertor.Lab
                {
                    alpha = bins[i].ac, L = bins[i].Lc, A = bins[i].Ac, B = bins[i].Bc
                };
                double alphaDiff = hasSemiTransparency ? Math.Abs(lab2.alpha - lab1.alpha) : 0;
                double nerr = nerr2 * Sqr(alphaDiff) * alphaDiff / 3.0;
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.L - lab1.L);
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.A - lab1.A);
                if (nerr >= err)
                    continue;

                nerr += (1 - ratio) * nerr2 * Sqr(lab2.B - lab1.B);

                if (nerr >= err)
                    continue;

                var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                nerr += ratio * nerr2 * Sqr(deltaL_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
                nerr += ratio * nerr2 * Sqr(deltaC_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
                nerr += ratio * nerr2 * Sqr(deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                nerr += ratio * nerr2 * CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
                if (nerr >= err)
                    continue;

                err = (float)nerr;
                nn = i;
            }
            bin1.err = (float) err;
            bin1.nn = nn;
        }
        protected override void Pnnquan(int[] pixels, Color[] palettes, int nMaxColors, bool quan_sqrt)
        {
            var bins = new Pnnbin[65536];

            /* Build histogram */
            foreach (var pixel in pixels)
            {
                // !!! Can throw gamma correction in here, but what to do about perceptual
                // !!! nonuniformity then?
                int index = GetARGBIndex(pixel, hasSemiTransparency);
                GetLab(pixel, out var lab1);
                if (bins[index] == null)
                    bins[index] = new Pnnbin();
                bins[index].ac += (float)lab1.alpha;
                bins[index].Lc += (float)lab1.L;
                bins[index].Ac += (float)lab1.A;
                bins[index].Bc += (float)lab1.B;
                bins[index].cnt++;
            }

            /* Cluster nonempty bins at one end of array */
            int maxbins = 0;
            var heap = new int[65537];
            int i = 0;
            for (; i < bins.Length; ++i)
            {
                if (bins[i] == null)
                    continue;

                var d = 1.0f / (float)bins[i].cnt;
                bins[i].ac *= d;
                bins[i].Lc *= d;
                bins[i].Ac *= d;
                bins[i].Bc *= d;

                bins[maxbins++] = bins[i];
            }

            var proportional = Sqr(nMaxColors) / maxbins;
            if ((proportional < .022 || proportional > .5) && nMaxColors < 64)
                quan_sqrt = false;
            
            for (i = 0; i < maxbins - 1; ++i)
            {
                bins[i].fw = i + 1;
                bins[i + 1].bk = i;

                if (quan_sqrt)
                    bins[i].cnt = (int) Math.Sqrt(bins[i].cnt);
            }
            if (quan_sqrt)
                bins[i].cnt = (int) Math.Sqrt(bins[i].cnt);

            int h, l, l2;
            ratio = 0.0;
            /* Initialize nearest neighbors and build heap of them */
            for (i = 0; i < maxbins; ++i)
            {
                Find_nn(bins, i);
                /* Push slot on heap */
                var err = bins[i].err;

                for (l = ++heap[0]; l > 1; l = l2)
                {
                    l2 = l >> 1;
                    if (bins[h = heap[l2]].err <= err)
                        break;
                    heap[l] = h;
                }
                heap[l] = i;
            }

            if (quan_sqrt && nMaxColors < 64)
                ratio = Math.Min(1.0, Math.Pow(nMaxColors, 1.55) / maxbins);
            else if (quan_sqrt)
                ratio = Math.Min(1.0, Math.Pow(nMaxColors, 1.05) / pixelMap.Count);            
            else
                ratio = .75;
            /* Merge bins which increase error the least */
            int extbins = maxbins - nMaxColors;
            for (i = 0; i < extbins;)
            {
                Pnnbin tb;
                /* Use heap to find which bins to merge */
                for (; ; )
                {
                    int b1 = heap[1];
                    tb = bins[b1]; /* One with least error */
                    /* Is stored error up to date? */
                    if ((tb.tm >= tb.mtm) && (bins[tb.nn].mtm <= tb.tm))
                        break;
                    if (tb.mtm == 0xFFFF) /* Deleted node */
                        b1 = heap[1] = heap[heap[0]--];
                    else /* Too old error value */
                    {
                        Find_nn(bins, b1);
                        tb.tm = i;
                    }
                    /* Push slot down */
                    var err = bins[b1].err;
                    for (l = 1; (l2 = l + l) <= heap[0]; l = l2)
                    {
                        if ((l2 < heap[0]) && (bins[heap[l2]].err > bins[heap[l2 + 1]].err))
                            ++l2;
                        if (err <= bins[h = heap[l2]].err)
                            break;
                        heap[l] = h;
                    }
                    heap[l] = b1;
                }

                /* Do a merge */
                var nb = bins[tb.nn];
                var n1 = tb.cnt;
                var n2 = nb.cnt;
                var d = 1.0f / (float)(n1 + n2);
                tb.ac = d * (n1 * tb.ac + n2 * nb.ac);
                tb.Lc = d * (n1 * tb.Lc + n2 * nb.Lc);
                tb.Ac = d * (n1 * tb.Ac + n2 * nb.Ac);
                tb.Bc = d * (n1 * tb.Bc + n2 * nb.Bc);
                tb.cnt += nb.cnt;
                tb.mtm = ++i;

                /* Unchain deleted bin */
                bins[nb.bk].fw = nb.fw;
                bins[nb.fw].bk = nb.bk;
                nb.mtm = 0xFFFF;
            }

            /* Fill palette */
            int k = 0;
            for (i = 0; ; ++k)
            {
                var lab1 = new CIELABConvertor.Lab
                {
                    alpha = bins[i].ac, L = bins[i].Lc, A = bins[i].Ac, B = bins[i].Bc
                };
                palettes[k] = CIELABConvertor.LAB2RGB(lab1);
                if (m_transparentPixelIndex >= 0 && palettes[k] == m_transparentColor)
                    Swap(ref palettes[0], ref palettes[k]);

                if ((i = bins[i].fw) == 0)
                    break;
            }
        }

        protected override ushort NearestColorIndex(Color[] palette, int nMaxColors, int argb)
        {
            if (nearestMap.TryGetValue(argb, out var k))
                return k;

            var c = Color.FromArgb(argb);

            double mindist = 1e100;
            GetLab(argb, out var lab1);

            for (int i = 0; i < nMaxColors; ++i)
            {
                var c2 = palette[i];

                var curdist = Sqr(c2.A - c.A);
                if (curdist > mindist)
                    continue;

                GetLab(c2.ToArgb(), out var lab2);
                if (nMaxColors > 32)
                {
                    curdist += PR * Sqr(c2.R - c.R);
                    if (curdist > mindist)
                        continue;

                    curdist += PG * Sqr(c2.G - c.G);
                    if (curdist > mindist)
                        continue;

                    curdist += PB * Sqr(c2.B - c.B);
                    if (PB < 1)
                    {
                        if (curdist > mindist)
                            continue;

                        curdist += Sqr(lab2.B - lab1.B) / 2.0;
                    }
                }
                else
                {  
                    var deltaL_prime_div_k_L_S_L = CIELABConvertor.L_prime_div_k_L_S_L(lab1, lab2);
                    curdist += Sqr(deltaL_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    var deltaC_prime_div_k_L_S_L = CIELABConvertor.C_prime_div_k_L_S_L(lab1, lab2, out var a1Prime, out var a2Prime, out var CPrime1, out var CPrime2);
                    curdist += Sqr(deltaC_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    var deltaH_prime_div_k_L_S_L = CIELABConvertor.H_prime_div_k_L_S_L(lab1, lab2, a1Prime, a2Prime, CPrime1, CPrime2, out var barCPrime, out var barhPrime);
                    curdist += Sqr(deltaH_prime_div_k_L_S_L);
                    if (curdist > mindist)
                        continue;

                    curdist += CIELABConvertor.R_T(barCPrime, barhPrime, deltaC_prime_div_k_L_S_L, deltaH_prime_div_k_L_S_L);
                }

                if (curdist > mindist)
                    continue;
                mindist = curdist;
                k = (ushort)i;
            }
            nearestMap[argb] = k;
            return k;
        }
        protected override ushort ClosestColorIndex(Color[] palette, int nMaxColors, int pixel)
        {
            ushort k = 0;
            if (!closestMap.TryGetValue(pixel, out var closest))
            {
                closest = new ushort[5];
                closest[2] = closest[3] = ushort.MaxValue;
                GetLab(pixel, out var lab1);

                for (; k < nMaxColors; ++k)
                {
                    var c2 = palette[k];
                    GetLab(c2.ToArgb(), out var lab2);

                    closest[4] = (ushort)(Sqr(lab2.alpha - lab1.alpha) + CIELABConvertor.CIEDE2000(lab2, lab1));
                    //closest[4] = (short) (Math.abs(lab2.alpha - lab1.alpha) + Math.abs(lab2.L - lab1.L) + Math.abs(lab2.A - lab1.A) + Math.abs(lab2.B - lab1.B));
                    if (closest[4] < closest[2])
                    {
                        closest[1] = closest[0];
                        closest[3] = closest[2];
                        closest[0] = (ushort)k;
                        closest[2] = closest[4];
                    }
                    else if (closest[4] < closest[3])
                    {
                        closest[1] = (ushort)k;
                        closest[3] = closest[4];
                    }
                }

                if (closest[3] == ushort.MaxValue)
                    closest[2] = 0;
            }

            if (closest[2] == 0 || (rand.Next(short.MaxValue) % (closest[3] + closest[2])) <= closest[3])
                k = closest[0];
            else
                k = closest[1];

            closestMap[pixel] = closest;
            return k;
        }       

    }
}