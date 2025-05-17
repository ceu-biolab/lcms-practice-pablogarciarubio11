package adduct;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.abs;

public class MassTransformation {

// !! TODO create functions to transform the mass of the mzs to monoisotopic masses and vice versa.
    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {
        Double monoMass=0.0;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).

        /*
        if Adduct is single charge the formula is M = m/z +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is M = ( mz +- adductMass ) * charge

        if adduct is a dimer or multimer the formula is M =  (mz +- adductMass) / numberOfMultimer

        return monoisotopicMass;

         */

        List<Double> infoAdduct=getInfoAdduct(adduct);
        Double multimer=infoAdduct.get(0);

        Double charge=abs(infoAdduct.get(1));
        monoMass = ((mz * charge) + infoAdduct.get(2)) / multimer;
        return monoMass;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoisotopicMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoisotopicMass, String adduct) {
        Double mz=0.0;
        // !! TODO METHOD
        // !! TODO Create the necessary regex to obtain the multimer (number before the M) and the charge (number before the + or - (if no number, the charge is 1).

        /*
        if Adduct is single charge the formula is m/z = M +- adductMass. Charge is 1 so it does not affect

        if Adduct is double or triple charged the formula is mz = M/charge +- adductMass

        if adduct is a dimer or multimer the formula is mz = M * numberOfMultimer +- adductMass

        return monoisotopicMass;

         */
        List<Double> infoAdduct=getInfoAdduct(adduct);
        Double multimer=infoAdduct.get(0);

        double charge=abs(infoAdduct.get(1));
        mz = ((monoisotopicMass * multimer) - infoAdduct.get(2))/charge;
        return mz;
    }

    //Method to get the information of the adduct:
    public static List<Double> getInfoAdduct(String adduct){
        List<Double> infoAdduct = new ArrayList<Double>();
        if (adduct == null) {return null;}
        /*
        The "-" sign character causes problems when comparing the string to the pattern, because there are 2 characters "-"
        that look the same but are different characters. So as a solution, we replace all the minus characters with one
        character minus.
         */
        String adduct1 = adduct.replace("−", "-").replace("–", "-").replace("—", "-");
        Pattern pattern = Pattern.compile("\\[(\\d*)M([+-].+)](\\d*[+-]?)");
        Matcher matcher = pattern.matcher(adduct1);

        if (!matcher.matches()) {
            throw new IllegalArgumentException("Formato de aducto inválido: " + adduct1);
        }

        int multimer = matcher.group(1).isEmpty() ? 1 : Integer.parseInt(matcher.group(1));
        infoAdduct.add(Double.valueOf(multimer));
        String chargeStr = matcher.group(3).trim();
        int charge = 1;
        if (chargeStr.matches("\\d+[+-]")) {
            charge = Integer.parseInt(chargeStr.substring(0, chargeStr.length() - 1));
            if (chargeStr.endsWith("-")) charge = -charge;
        } else if (chargeStr.equals("-")) {
            charge = -1;
        }
        infoAdduct.add(Double.valueOf(charge));
        Map<String, Double> searchList = new LinkedHashMap<>();
        if(infoAdduct.get(1)>0){
            searchList.putAll(AdductList.MAPMZPOSITIVEADDUCTS);
        }else{
            searchList.putAll(AdductList.MAPMZNEGATIVEADDUCTS);
        }

        /*
        Due to the replacement of the minus signs that I have done before, when I search for the adduct mass,
        the keys contain the minus sign that I have replaced. To solve this, I convert the string back to its
        original form.
         */
        // String originalAdduct = adduct.replace("-", "−");
        if(searchList.containsKey(adduct)){
            infoAdduct.add(searchList.get(adduct));
        }else{
            infoAdduct.add(0.0);
        }

        return infoAdduct;
    }


    public static void main(String[] args) {
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[M+H]+"));
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[2M+H]+"));
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[M+2H]2+"));
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[M+Cl]−"));
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[M-H]−"));
        System.out.println(getMonoisotopicMassFromMZ(565.3,"[M-2H]2−"));


        System.out.println("");

        System.out.println(getMZFromMonoisotopicMass(564.2927239999999,"[M+H]+"));
        System.out.println(getMZFromMonoisotopicMass(282.14636199999995,"[2M+H]+"));
        System.out.println(getMZFromMonoisotopicMass(1128.5854479999998,"[M+2H]2+"));
        System.out.println(getMZFromMonoisotopicMass(530.330598,"[M+Cl]−"));
        System.out.println(getMZFromMonoisotopicMass(566.307276,"[M-H]−"));
        System.out.println(getMZFromMonoisotopicMass(1131.607276,"[M-2H]2−"));

    }
}

