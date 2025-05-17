package lipid;

import adduct.AdductList;
import adduct.MassTransformation;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IoniationMode ionizationMode;
    private String adduct; // !!TODO The adduct will be detected based on the groupedSignals
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IoniationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        // !!TODO This set should be sorted according to help the program to deisotope the signals plus detect the adduct
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        this.adduct = this.detectAdduct();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IoniationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    // !!TODO Detect the adduct with an algorithm or with drools, up to the user.



    public String detectAdduct() { //FINAL FUNCTION
        double mzTolerance = 0.2;
        if (groupedSignals == null || groupedSignals.size() < 2) {
            System.out.println("detectAdduct: Not enough peaks (" +
                    (groupedSignals == null ? 0 : groupedSignals.size()) + ")");
            return "Unknown";
        }

        Map<String, Double> adductMap = new LinkedHashMap<>();
        if(getIonizationMode()==IoniationMode.POSITIVE){
            adductMap=AdductList.MAPMZPOSITIVEADDUCTS;
        }else{
            adductMap=AdductList.MAPMZNEGATIVEADDUCTS;
        }

        double observedMz = this.getMz();
        System.out.println("detectAdduct: observedMz = " + observedMz + ", mode = " + getIonizationMode());

        //For each candidateAdduct that could explain observeMz:
        for (String candidateAdduct : adductMap.keySet()) {
            System.out.println("  Testing candidateAdduct: " + candidateAdduct);
            try {
                //I calculate the monoisotopic mass according to the candidateAdduct for observedMz
                double monoisotopicMass = MassTransformation.getMonoisotopicMassFromMZ(observedMz, candidateAdduct);
                System.out.println("    monoisotopicMass = " + monoisotopicMass);

                //I search so that some of the other peaks corresponds to that mass
                for (Peak otherPeak : groupedSignals) {
                    System.out.println("    Comparing with otherPeak: " + otherPeak);
                    if (Math.abs(otherPeak.getMz() - observedMz) <= mzTolerance) {
                        //It's the same peak as the objective, so I skip it.
                        System.out.println("      Skip: same as observedMz");
                        continue;
                    }
                    //I try all the possible adducts for this otherPeak
                    for (String secondAdduct : adductMap.keySet()) {
                        double expectedMz = MassTransformation
                                .getMZFromMonoisotopicMass(monoisotopicMass, secondAdduct);
                        double diff = Math.abs(expectedMz - otherPeak.getMz());
                        System.out.println("      secondAdduct=" + secondAdduct + ", expectedMz=" + expectedMz + ", observed=" + otherPeak.getMz() + ", diff=" + diff);
                        if (diff <= mzTolerance) {
                            System.out.println("    DETECTED adduct: " + candidateAdduct + " (via " + secondAdduct + ")");
                            return candidateAdduct;
                        }
                    }
                }

            } catch (IllegalArgumentException e) {
                //To catch some error in the parsing of the adduct
                System.out.println("    Ignored candidateAdduct (parse error): " + candidateAdduct);
            }
        }

        System.out.println("detectAdduct: No adduct detected");
        return "Unknown";

    }



    public static void main(String[] args) {
        Peak mH = new Peak(700.500, 80000.0); // [M+H]+
        Peak mNa = new Peak(722.482, 100000.0);  // [M+Na]+ DETECT THIS ONE
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);

        double annotationMZ = 722.482d;
        double annotationIntensity = 100000.0;
        double annotationRT = 6.5d;
        Annotation annotation = new Annotation(lipid, annotationMZ, annotationIntensity, annotationRT, IoniationMode.POSITIVE, Set.of(mH, mNa));
        System.out.println(annotation);

        System.out.println("___________________________________");

        Peak singlyCharged = new Peak(700.500, 100000.0);  // [M+H]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+ DETECT THIS ONE

        Lipid lipid1 = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
        Annotation annotation1 = new Annotation(lipid1, doublyCharged.getMz(), doublyCharged.getIntensity(), 10d, IoniationMode.POSITIVE, Set.of(singlyCharged, doublyCharged));
        System.out.println(annotation1);

        System.out.println("___________________________________");

        Peak mH2 = new Peak(700.500, 100000.0); // [M+H]+
        Peak mk2 = new Peak(738.4564, 80000.0);  // [M+K]+ DETECT THIS ONE
        Lipid lipid2 = new Lipid(1, "PC 34:1", "C42H82NO8P", "PC", 34, 1);

        double annotationRT2 = 6.5d;
        Annotation annotation2 = new Annotation(lipid, mk2.getMz(), mk2.getIntensity(), annotationRT2, IoniationMode.POSITIVE, Set.of(mH2, mk2));
        System.out.println(annotation2);

        System.out.println("___________________________________");

        Peak peak1 = new Peak(350.754, 85000.0);   // [M+2H]2+
        Peak peak2 = new Peak(1399.992724, 85000.0);   // [2M+H]+ DETECT THIS ONE

        Lipid lipid3 = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
        Annotation annotation3 = new Annotation(lipid3, peak2.getMz(), peak2.getIntensity(), 10d, IoniationMode.POSITIVE, Set.of(peak1,peak2));
        System.out.println(annotation3);

        System.out.println("___________________________________");

        Peak mh1    = new Peak(883.7760,  85000.0);  // [M–H]⁻ DETECT THIS ONE
        Peak mCl    = new Peak(919.7522,  85000.0);  // [M+Cl]⁻
        Lipid lipid4 = new Lipid(3, "TG 54:3", "C57H104O6", "TG", 54, 3);
        Annotation annotation4 = new Annotation(lipid4,mh1.getMz(),mh1.getIntensity(),10d,IoniationMode.NEGATIVE,Set.of(mh1, mCl));
    }
}
