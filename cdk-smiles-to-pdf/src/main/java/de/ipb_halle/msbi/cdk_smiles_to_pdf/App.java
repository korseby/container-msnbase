package de.ipb_halle.msbi.cdk_smiles_to_pdf;

import java.io.IOException;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Hello world!
 *
 */
public class App {
	public static void main(String[] args) throws IOException, CDKException {

		if (args.length == 2) {

			String smiles = args[0];
			String outfile = args[1];

			IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
			SmilesParser smipar = new SmilesParser(bldr);

			IAtomContainer mol = smipar.parseSmiles(smiles);
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
