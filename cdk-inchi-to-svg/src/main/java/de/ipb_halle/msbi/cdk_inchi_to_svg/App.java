package de.ipb_halle.msbi.cdk_inchi_to_svg;

import java.io.IOException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import net.sf.jniinchi.INCHI_RET;
/**
 * Hello world!
 *
 */
public class InChitoSVG {
	public static void main(String[] args) throws IOException, CDKException {
		if (args.length == 2) {
			String inchi = args[0];
			String outfile = args[1];
			
			
			// Get InChIToStructure
			InChIToStructure intostruct = InChIGeneratorFactory.getInstance().getInChIToStructure(inchi, DefaultChemObjectBuilder.getInstance());
			INCHI_RET ret = intostruct.getReturnStatus();
			if (ret == INCHI_RET.WARNING) {
				// Structure generated, but with warning message
				System.out.println("InChI warning: " + intostruct.getMessage());
			}
			IAtomContainer mol = intostruct.getAtomContainer();
			//mol.setProperty(CDKConstants.TITLE, "caffeine"); // title already set from input!
			DepictionGenerator dptgen = new DepictionGenerator();
			dptgen.withAtomColors()
					// .withSize(200, 250) // px (raster) or mm (vector)
					// .withMolTitle()
					// .withTitleColor(Color.DARK_GRAY) // annotations are red by default
					.depict(mol).writeTo(outfile); // svg preferred but will be 20cmx25cm!
		}
	}
}
