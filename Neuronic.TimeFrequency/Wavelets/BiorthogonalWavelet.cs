using Accord.Math;
using System;

namespace Neuronic.TimeFrequency.Wavelets
{
    public class BiorthogonalWavelet : OrthogonalWavelet
    {
        public BiorthogonalWavelet(string shortName, string familyName, double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, double freq = 0)
             : base(shortName, familyName, lowRec, highRec, lowDec, highDec, vanishingMoments, freq)
        {
            Other = new BiorthogonalWavelet(this, shortName, familyName, lowDec.Reversed(), highDec.Reversed(), lowRec.Reversed(), highRec.Reversed(), vanishingMoments, freq);
        }

        protected BiorthogonalWavelet(BiorthogonalWavelet reverse, string shortName, string familyName, double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, double freq = 0)
             : base(shortName, familyName, lowRec, highRec, lowDec, highDec, vanishingMoments, freq)
        {
            Other = reverse;
        }

        public BiorthogonalWavelet Other { get; private set; }

        protected override int DesiredOutputSize => base.DesiredOutputSize - 1;

        protected override double[] Upcoef(double[] coeffs, int level)
        {
            if (VanishingMoments % 4 != 1)
                coeffs[0] *= -1;
            return base.Upcoef(coeffs, level);
        }
    }
}
