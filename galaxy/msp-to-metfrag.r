<?xml version='1.0' encoding='UTF-8'?>
<tool id="msp-to-metfrag" name="msp-to-metfrag" version="0.1">
  <requirements>
    <container type="docker">korseby/mtbls709</container>
  </requirements>
  <description>Read MS2 spectra from MSP file and create MetFrag parameter file.</description>
  <stdio>
    <regex match="" source="stderr" level="warning" description="R messages" />
    <exit_code range="1:" level="fatal" description="Tool error" />
  </stdio>
  <command><![CDATA[
/usr/local/bin/msp-to-metfrag.r $mspfile $metfragparameterfile;
  ]]>
  </command>
  <inputs>
    <param name="mspfile" type="data" format="txt" optional="False" label="MSP file" help="MSP file" />
  </inputs>
  <outputs>
    <data name="metfragparameterfile" type="data" format="txt" label="${mspfile.display_name}_parameter.txt" />
  </outputs>
  <help>
.. class:: infomark

**Authors**

| **Kristian Peters (kpeters@ipb-halle.de)** wrote and maintains this module.

---------------------------------------------------

========================
MSP to MetFrag parameter
========================

-----------
Description
-----------

| Read spectra from MSP file and create MetFrag parameter file.

-----------
Input files
-----------

+------------------------------+------------+
| File                         |   Format   |
+==============================+============+
| 1)  MSP file(s)              |   MSP      |
+------------------------------+------------+

----------
Parameters
----------
	  
XCMS-Set file
	| A rdata file with XCMS-Set objects that were grouped by **xcms-collect-peaks**
        |

------------
Output files
------------
	
xcmssplit_input_1.rdata
	| (Multiple) rdata file(s) containing one XCMS-Set object each
        |


---------------------------------------------------

  </help>
</tool>

