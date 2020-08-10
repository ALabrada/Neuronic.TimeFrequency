namespace Neuronic.TimeFrequency.Wavelets
{
    public class BiorthogonalWavelet : OrthogonalWavelet
    {
        public BiorthogonalWavelet(double[] lowRec, double[] highRec, double[] lowDec, double[] highDec, int vanishingMoments, OrthogonalWavelet other, double freq = 0) : base(lowRec, highRec, lowDec, highDec, vanishingMoments, freq)
        {
        }

        public OrthogonalWavelet Other { get; }

        protected override double[] Upcoef(double[] coeffs, int level, bool recA = false)
        {
            if (VanishingMoments % 4 != 1)
                coeffs[0] *= -1;
            return base.Upcoef(coeffs, level, recA);
        }
    }
}
