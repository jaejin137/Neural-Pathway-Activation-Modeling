function out = model
%
% Barb_R_STNGPeDBS_rev1_ver3_STN-1-3-0p45_solved.m
%
% Model exported on Feb 14 2024, 15:40 by COMSOL 6.1.0.357.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\COMSOL_Output');

model.label('Barb_R_STNGPeDBS_rev1_ver3_STN-1-3-0p45_solved.mph');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component.create('mcomp1', 'MeshComponent');

model.geom.create('mgeom1', 3);

model.component.create('mcomp2', 'MeshComponent');

model.geom.create('mgeom2', 3);

model.component.create('mcomp3', 'MeshComponent');

model.geom.create('mgeom3', 3);

model.component.create('mcomp4', 'MeshComponent');

model.geom.create('mgeom4', 3);

model.component.create('mcomp5', 'MeshComponent');

model.geom.create('mgeom5', 3);

model.component.create('mcomp6', 'MeshComponent');

model.geom.create('mgeom6', 3);

model.component.create('mcomp7', 'MeshComponent');

model.geom.create('mgeom7', 3);

model.component.create('mcomp8', 'MeshComponent');

model.geom.create('mgeom8', 3);

model.component.create('mcomp9', 'MeshComponent');

model.geom.create('mgeom9', 3);

model.component.create('mcomp10', 'MeshComponent');

model.geom.create('mgeom10', 3);

model.component.create('mcomp11', 'MeshComponent');

model.geom.create('mgeom11', 3);

model.component.create('mcomp12', 'MeshComponent');

model.geom.create('mgeom12', 3);

model.component.create('mcomp13', 'MeshComponent');

model.geom.create('mgeom13', 3);

model.component.create('mcomp14', 'MeshComponent');

model.geom.create('mgeom14', 3);

model.component.create('mcomp15', 'MeshComponent');

model.geom.create('mgeom15', 3);

model.component.create('mcomp16', 'MeshComponent');

model.geom.create('mgeom16', 3);

model.component.create('mcomp17', 'MeshComponent');

model.geom.create('mgeom17', 3);

model.component.create('mcomp18', 'MeshComponent');

model.geom.create('mgeom18', 3);

model.component.create('mcomp19', 'MeshComponent');

model.geom.create('mgeom19', 3);

model.component.create('mcomp20', 'MeshComponent');

model.geom.create('mgeom20', 3);

model.component.create('mcomp21', 'MeshComponent');

model.geom.create('mgeom21', 3);

model.component.create('mcomp22', 'MeshComponent');

model.geom.create('mgeom22', 3);

model.component.create('mcomp23', 'MeshComponent');

model.geom.create('mgeom23', 3);

model.component.create('mcomp24', 'MeshComponent');

model.geom.create('mgeom24', 3);

model.component.create('mcomp25', 'MeshComponent');

model.geom.create('mgeom25', 3);

model.component('comp1').curvedInterior(false);

model.result.table.create('evl3', 'Table');

model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int1').set('source', 'file');
model.component('comp1').func('int1').label('conductivity');
model.component('comp1').func('int1').set('importedname', 'nv_3049.txt');
model.component('comp1').func('int1').set('importedstruct', 'Spreadsheet');
model.component('comp1').func('int1').set('importeddim', '3D');
model.component('comp1').func('int1').set('funcs', {'xx' '1';  ...
'xy' '2';  ...
'xz' '3';  ...
'yy' '4';  ...
'yz' '5';  ...
'zz' '6'});
model.component('comp1').func('int1').set('defvars', true);
model.component('comp1').func('int1').set('frame', 'mesh');
model.component('comp1').func('int1').set('interp', 'neighbor');
model.component('comp1').func('int1').set('extrap', 'value');
model.component('comp1').func('int1').set('extrapvalue', '1e-6');
model.component('comp1').func('int1').set('fununit', {'S/m' 'S/m' 'S/m' 'S/m' 'S/m' 'S/m'});
model.component('comp1').func('int1').set('argunit', {'mm' 'mm' 'mm'});
model.component('comp1').func('int1').set('filename', 'C:\Users\jaejin\AppData\Local\Temp\2\csmphserver46534\nv_3049.txt');
model.component('comp1').func('int1').importData;
model.component('comp1').func('int1').set('nargs', '3');
model.component('comp1').func('int1').set('struct', 'spreadsheet');
model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').label('permittivity');
model.component('comp1').func('int2').set('importedname', 'nv_3049.txt');
model.component('comp1').func('int2').set('importedstruct', 'Spreadsheet');
model.component('comp1').func('int2').set('importeddim', '3D');
model.component('comp1').func('int2').set('funcs', {'ee' '7'});
model.component('comp1').func('int2').set('defvars', true);
model.component('comp1').func('int2').set('frame', 'mesh');
model.component('comp1').func('int2').set('interp', 'neighbor');
model.component('comp1').func('int2').set('extrap', 'value');
model.component('comp1').func('int2').set('extrapvalue', '1e-6');
model.component('comp1').func('int2').set('fununit', {'1'});
model.component('comp1').func('int2').set('argunit', {'mm' 'mm' 'mm'});
model.component('comp1').func('int2').set('filename', 'C:\Users\jaejin\AppData\Local\Temp\2\csmphserver46534\18b2s4e01thbz\nv_3049.txt');
model.component('comp1').func('int2').importData;
model.component('comp1').func('int2').set('nargs', '3');
model.component('comp1').func('int2').set('struct', 'spreadsheet');

model.component('comp1').mesh.create('mesh1');
model.mesh.create('mpart1', 'mgeom1');
model.mesh.create('mpart2', 'mgeom2');
model.mesh.create('mpart3', 'mgeom3');
model.mesh.create('mpart4', 'mgeom4');
model.mesh.create('mpart5', 'mgeom5');
model.mesh.create('mpart6', 'mgeom6');
model.mesh.create('mpart7', 'mgeom7');
model.mesh.create('mpart8', 'mgeom8');
model.mesh.create('mpart9', 'mgeom9');
model.mesh.create('mpart10', 'mgeom10');
model.mesh.create('mpart11', 'mgeom11');
model.mesh.create('mpart12', 'mgeom12');
model.mesh.create('mpart13', 'mgeom13');
model.mesh.create('mpart14', 'mgeom14');
model.mesh.create('mpart15', 'mgeom15');
model.mesh.create('mpart16', 'mgeom16');
model.mesh.create('mpart17', 'mgeom17');
model.mesh.create('mpart18', 'mgeom18');
model.mesh.create('mpart19', 'mgeom19');
model.mesh.create('mpart20', 'mgeom20');
model.mesh.create('mpart21', 'mgeom21');
model.mesh.create('mpart22', 'mgeom22');
model.mesh.create('mpart23', 'mgeom23');
model.mesh.create('mpart24', 'mgeom24');
model.mesh.create('mpart25', 'mgeom25');

model.geom('mgeom1').lengthUnit('mm');

model.mesh('mpart1').create('imp1', 'Import');
model.mesh('mpart1').feature('imp1').set('source', 'stl');
model.mesh('mpart1').feature('imp1').set('filename', 'C:\Users\slops001\Box\DBS_Patient_Modeling\Essential_Tremor\ET025\DTI_brain_mask7.stl');
model.mesh('mpart1').feature('imp1').importData;
model.mesh('mpart1').run;

model.geom('mgeom2').lengthUnit('mm');

model.mesh('mpart2').create('imp1', 'Import');
model.mesh('mpart2').feature('imp1').set('source', 'stl');
model.mesh('mpart2').feature('imp1').set('filename', 'C:\Users\slops001\Box\DBS_Patient_Modeling\COMSOL templates\Skull_little.stl');
model.mesh('mpart2').feature('imp1').importData;
model.mesh('mpart2').run;

model.geom('mgeom3').lengthUnit('mm');

model.mesh('mpart3').create('imp1', 'Import');
model.mesh('mpart3').feature('imp1').set('source', 'stl');
model.mesh('mpart3').feature('imp1').set('filename', 'C:\Users\kovac149\Box\DBS_Patient_Modeling\Essential_Tremor\ET024\ET024_brain_mask_DTI.stl');
model.mesh('mpart3').feature('imp1').importData;
model.mesh('mpart3').run;

model.geom('mgeom4').lengthUnit('mm');

model.mesh('mpart4').create('imp1', 'Import');
model.mesh('mpart4').feature('imp1').set('source', 'stl');
model.mesh('mpart4').feature('imp1').set('filename', 'C:\Users\kovac149\Box\DBS_Patient_Modeling\Essential_Tremor\ET024\ET024_brain_mask_DTI_clean.stl');
model.mesh('mpart4').feature('imp1').importData;
model.mesh('mpart4').run;

model.geom('mgeom5').lengthUnit('mm');

model.mesh('mpart5').create('imp1', 'Import');
model.mesh('mpart5').feature('imp1').set('source', 'stl');
model.mesh('mpart5').feature('imp1').set('filename', 'C:\Users\kovac149\Box\DBS_Patient_Modeling\Essential_Tremor\ET024\ET024_brain_mask_DTI_clean.stl');
model.mesh('mpart5').feature('imp1').importData;
model.mesh('mpart5').run;

model.geom('mgeom6').lengthUnit('mm');

model.mesh('mpart6').create('imp1', 'Import');
model.mesh('mpart6').feature('imp1').set('source', 'stl');
model.mesh('mpart6').feature('imp1').set('filename', 'C:\Users\kovac149\Box\DBS_Patient_Modeling\Essential_Tremor\ET023\DTI_brain_mask.stl');
model.mesh('mpart6').feature('imp1').importData;
model.mesh('mpart6').run;

model.geom('mgeom7').lengthUnit('mm');

model.mesh('mpart7').create('imp1', 'Import');
model.mesh('mpart7').feature('imp1').set('source', 'stl');
model.mesh('mpart7').feature('imp1').set('filename', 'C:\Users\chuxx366\Box\DBS_Patient_Modeling\Essential_Tremor\ET023\Brain_DTI.stl');
model.mesh('mpart7').feature('imp1').importData;
model.mesh('mpart7').run;

model.geom('mgeom8').lengthUnit('mm');

model.mesh('mpart8').create('imp1', 'Import');
model.mesh('mpart8').feature('imp1').set('source', 'stl');
model.mesh('mpart8').feature('imp1').set('filename', 'C:\Users\chuxx366\Box\DBS_Patient_Modeling\Essential_Tremor\ET023\Brain_DTI_1.stl');
model.mesh('mpart8').feature('imp1').importData;
model.mesh('mpart8').run;

model.geom('mgeom9').lengthUnit('mm');

model.mesh('mpart9').create('imp1', 'Import');
model.mesh('mpart9').feature('imp1').set('source', 'stl');
model.mesh('mpart9').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\UD1015_PD088\Data Processing\UD1015_BET_mask_DTIspace_expanded_smooth.stl');
model.mesh('mpart9').feature('imp1').importData;
model.mesh('mpart9').run;

model.geom('mgeom10').lengthUnit('mm');

model.mesh('mpart10').create('imp1', 'Import');
model.mesh('mpart10').feature('imp1').set('source', 'stl');
model.mesh('mpart10').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\Skull_little.stl');
model.mesh('mpart10').feature('imp1').importData;
model.mesh('mpart10').run;

model.geom('mgeom11').lengthUnit('mm');

model.mesh('mpart11').create('imp1', 'Import');
model.mesh('mpart11').feature('imp1').set('source', 'stl');
model.mesh('mpart11').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\UD1015_PD088\Data Processing\UD1015_BET_mask_DTIspace_expanded_smooth.stl');
model.mesh('mpart11').feature('imp1').importData;
model.mesh('mpart11').run;

model.geom('mgeom12').lengthUnit('mm');

model.mesh('mpart12').create('imp1', 'Import');
model.mesh('mpart12').feature('imp1').set('source', 'stl');
model.mesh('mpart12').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\COMSOL_Files\Skull_little.stl');
model.mesh('mpart12').feature('imp1').importData;
model.mesh('mpart12').run;

model.geom('mgeom13').lengthUnit('mm');

model.mesh('mpart13').create('imp1', 'Import');
model.mesh('mpart13').feature('imp1').set('source', 'stl');
model.mesh('mpart13').feature('imp1').set('filename', 'C:\Users\linnx092\Box\COMSOL templates\Medtronic_3389\MDT_3389_lead_200mm_tall_0p15mm_encapsulation_only_08_21_2017 (short).stl');
model.mesh('mpart13').feature('imp1').importData;
model.mesh('mpart13').run;

model.geom('mgeom14').lengthUnit('mm');

model.mesh('mpart14').create('imp1', 'Import');
model.mesh('mpart14').feature('imp1').set('source', 'stl');
model.mesh('mpart14').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\UD1001_CFR001\Data Processing\COMSOL Input\UD1001_BET_mask_DTIspace_expanded_smooth.stl');
model.mesh('mpart14').feature('imp1').importData;
model.mesh('mpart14').run;

model.geom('mgeom15').lengthUnit('mm');

model.mesh('mpart15').create('imp1', 'Import');
model.mesh('mpart15').feature('imp1').set('source', 'stl');
model.mesh('mpart15').feature('imp1').set('filename', 'C:\Users\linnx092\Box\Globus_Pallidus\COMSOL_Files\Skull_little.stl');
model.mesh('mpart15').feature('imp1').importData;
model.mesh('mpart15').run;

model.geom('mgeom16').lengthUnit('mm');

model.mesh('mpart16').create('imp1', 'Import');
model.mesh('mpart16').feature('imp1').set('source', 'stl');
model.mesh('mpart16').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data\Segmentations\STN_R.stl');
model.mesh('mpart16').feature('imp1').set('createdom', true);
model.mesh('mpart16').feature('imp1').importData;
model.mesh('mpart16').run;

model.geom('mgeom17').lengthUnit('mm');

model.mesh('mpart17').create('imp1', 'Import');
model.mesh('mpart17').feature('imp1').set('source', 'stl');
model.mesh('mpart17').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data\Segmentations\GPe_R.stl');
model.mesh('mpart17').feature('imp1').set('createdom', true);
model.mesh('mpart17').feature('imp1').importData;
model.mesh('mpart17').run;

model.geom('mgeom18').lengthUnit('mm');

model.mesh('mpart18').create('imp1', 'Import');
model.mesh('mpart18').feature('imp1').set('source', 'stl');
model.mesh('mpart18').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data\Segmentations\GPi_R.stl');
model.mesh('mpart18').feature('imp1').set('createdom', true);
model.mesh('mpart18').feature('imp1').importData;
model.mesh('mpart18').run;

model.geom('mgeom19').lengthUnit('mm');

model.mesh('mpart19').create('imp1', 'Import');
model.mesh('mpart19').feature('imp1').set('source', 'stl');
model.mesh('mpart19').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data\STL Files\brain_whole.stl');
model.mesh('mpart19').feature('imp1').set('createdom', true);
model.mesh('mpart19').feature('imp1').importData;
model.mesh('mpart19').run;

model.geom('mgeom20').lengthUnit('mm');

model.mesh('mpart20').create('imp1', 'Import');
model.mesh('mpart20').feature('imp1').set('source', 'stl');
model.mesh('mpart20').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\STL\STN_R.stl');
model.mesh('mpart20').feature('imp1').set('createdom', true);
model.mesh('mpart20').feature('imp1').importData;
model.mesh('mpart20').run;

model.geom('mgeom21').lengthUnit('mm');

model.mesh('mpart21').create('imp1', 'Import');
model.mesh('mpart21').feature('imp1').set('source', 'stl');
model.mesh('mpart21').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\STL\GPe_R.stl');
model.mesh('mpart21').feature('imp1').set('createdom', true);
model.mesh('mpart21').feature('imp1').importData;
model.mesh('mpart21').run;

model.geom('mgeom22').lengthUnit('mm');

model.mesh('mpart22').create('imp1', 'Import');
model.mesh('mpart22').feature('imp1').set('source', 'stl');
model.mesh('mpart22').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Segmentations\GPe_R.stl');
model.mesh('mpart22').feature('imp1').set('createdom', true);
model.mesh('mpart22').feature('imp1').importData;
model.mesh('mpart22').run;

model.geom('mgeom23').lengthUnit('mm');

model.mesh('mpart23').create('imp1', 'Import');
model.mesh('mpart23').feature('imp1').set('source', 'stl');
model.mesh('mpart23').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data\STL Files\GPi_R.stl');
model.mesh('mpart23').feature('imp1').set('createdom', true);
model.mesh('mpart23').feature('imp1').importData;
model.mesh('mpart23').run;

model.geom('mgeom24').lengthUnit('mm');

model.mesh('mpart24').create('imp1', 'Import');
model.mesh('mpart24').feature('imp1').set('source', 'stl');
model.mesh('mpart24').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\STL\Cerebrum_smoothed.stl');
model.mesh('mpart24').feature('imp1').set('createdom', true);
model.mesh('mpart24').feature('imp1').importData;
model.mesh('mpart24').run;

model.geom('mgeom25').lengthUnit('mm');

model.mesh('mpart25').create('imp1', 'Import');
model.mesh('mpart25').feature('imp1').set('source', 'stl');
model.mesh('mpart25').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Segmentations\brain_whole.stl');
model.mesh('mpart25').feature('imp1').set('createdom', true);
model.mesh('mpart25').feature('imp1').importData;
model.mesh('mpart25').run;

model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').geomRep('comsol');
model.component('comp1').geom('geom1').create('imp8', 'Import');
model.component('comp1').geom('geom1').feature('imp8').label('Import_STN_Right');
model.component('comp1').geom('geom1').feature('imp8').set('type', 'mesh');
model.component('comp1').geom('geom1').feature('imp8').set('mesh', 'mpart20');
model.component('comp1').geom('geom1').feature('imp8').importData;
model.component('comp1').geom('geom1').create('imp9', 'Import');
model.component('comp1').geom('geom1').feature('imp9').label('Import_GPe_Right');
model.component('comp1').geom('geom1').feature('imp9').set('type', 'mesh');
model.component('comp1').geom('geom1').feature('imp9').set('mesh', 'mpart22');
model.component('comp1').geom('geom1').feature('imp9').importData;
model.component('comp1').geom('geom1').create('imp10', 'Import');
model.component('comp1').geom('geom1').feature('imp10').label('Import_GPi_Right');
model.component('comp1').geom('geom1').feature('imp10').set('type', 'mesh');
model.component('comp1').geom('geom1').feature('imp10').set('mesh', 'mpart23');
model.component('comp1').geom('geom1').feature('imp10').importData;
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').label('Import_Heraeus_lead_STN_R');
model.component('comp1').geom('geom1').feature('imp1').set('type', 'native');
model.component('comp1').geom('geom1').feature('imp1').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\COMSOL_Files_JL\Heraeus_40mm_0p25mm_encap_minus_lead.mphbin');
model.component('comp1').geom('geom1').create('imp2', 'Import');
model.component('comp1').geom('geom1').feature('imp2').label('Import_encapsulation');
model.component('comp1').geom('geom1').feature('imp2').set('type', 'native');
model.component('comp1').geom('geom1').feature('imp2').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\COMSOL_Files_JL\Heraeus_40mm_0p25mm_encap_only.mphbin');
model.component('comp1').geom('geom1').create('rot6', 'Rotate');
model.component('comp1').geom('geom1').feature('rot6').set('specify', 'eulerang');
model.component('comp1').geom('geom1').feature('rot6').set('eulerang', [30.4 -35.7 40]);
model.component('comp1').geom('geom1').feature('rot6').selection('input').set({'imp1' 'imp2'});
model.component('comp1').geom('geom1').create('mov8', 'Move');
model.component('comp1').geom('geom1').feature('mov8').setIndex('displx', '-6.78', 0);
model.component('comp1').geom('geom1').feature('mov8').setIndex('disply', '1.15', 0);
model.component('comp1').geom('geom1').feature('mov8').setIndex('displz', '23', 0);
model.component('comp1').geom('geom1').feature('mov8').selection('input').set({'rot6'});
model.component('comp1').geom('geom1').create('imp11', 'Import');
model.component('comp1').geom('geom1').feature('imp11').label('Import_Heraeus_lead_GPe_R');
model.component('comp1').geom('geom1').feature('imp11').set('type', 'native');
model.component('comp1').geom('geom1').feature('imp11').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\COMSOL_Files_JL\Heraeus_40mm_0p25mm_encap_minus_lead.mphbin');
model.component('comp1').geom('geom1').create('imp12', 'Import');
model.component('comp1').geom('geom1').feature('imp12').label('Import_encapsulation 1');
model.component('comp1').geom('geom1').feature('imp12').set('type', 'native');
model.component('comp1').geom('geom1').feature('imp12').set('filename', 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\COMSOL_Files_JL\Heraeus_40mm_0p25mm_encap_only.mphbin');
model.component('comp1').geom('geom1').create('rot7', 'Rotate');
model.component('comp1').geom('geom1').feature('rot7').set('specify', 'eulerang');
model.component('comp1').geom('geom1').feature('rot7').set('eulerang', [89 -48.8 40]);
model.component('comp1').geom('geom1').feature('rot7').selection('input').set({'imp11' 'imp12'});
model.component('comp1').geom('geom1').create('mov9', 'Move');
model.component('comp1').geom('geom1').feature('mov9').setIndex('displx', '-11.29', 0);
model.component('comp1').geom('geom1').feature('mov9').setIndex('disply', '-0.8', 0);
model.component('comp1').geom('geom1').feature('mov9').setIndex('displz', '28.54', 0);
model.component('comp1').geom('geom1').feature('mov9').selection('input').set({'rot7'});
model.component('comp1').geom('geom1').create('imp7', 'Import');
model.component('comp1').geom('geom1').feature('imp7').label('Brain_Import');
model.component('comp1').geom('geom1').feature('imp7').set('type', 'mesh');
model.component('comp1').geom('geom1').feature('imp7').set('mesh', 'mpart25');
model.component('comp1').geom('geom1').feature('imp7').importData;
model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
model.component('comp1').geom('geom1').feature('cyl1').set('pos', [0 0 -30]);
model.component('comp1').geom('geom1').feature('cyl1').set('r', 50);
model.component('comp1').geom('geom1').feature('cyl1').set('h', 100);
model.component('comp1').geom('geom1').create('rot3', 'Rotate');
model.component('comp1').geom('geom1').feature('rot3').setIndex('rot', '45', 0);
model.component('comp1').geom('geom1').feature('rot3').selection('input').set({'cyl1'});
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
model.component('comp1').geom('geom1').feature('sph1').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('sph1').set('r', 100);
model.component('comp1').geom('geom1').create('dif3', 'Difference');
model.component('comp1').geom('geom1').feature('dif3').label('Difference 3 (ext_sph - head) ');
model.component('comp1').geom('geom1').feature('dif3').set('keepadd', true);
model.component('comp1').geom('geom1').feature('dif3').set('keepsubtract', true);
model.component('comp1').geom('geom1').feature('dif3').selection('input').set({'sph1(1)'});
model.component('comp1').geom('geom1').feature('dif3').selection('input2').set({'rot3(1)'});
model.component('comp1').geom('geom1').create('dif2', 'Difference');
model.component('comp1').geom('geom1').feature('dif2').label('Difference 2 (encap/lead - ext_sph)');
model.component('comp1').geom('geom1').feature('dif2').selection('input').set({'mov8(1)' 'mov9(1)'});
model.component('comp1').geom('geom1').feature('dif2').selection('input2').set({'dif3'});
model.component('comp1').geom('geom1').create('dif8', 'Difference');
model.component('comp1').geom('geom1').feature('dif8').label('Difference 8(head - brain/encap)');
model.component('comp1').geom('geom1').feature('dif8').set('keepadd', true);
model.component('comp1').geom('geom1').feature('dif8').set('keepsubtract', true);
model.component('comp1').geom('geom1').feature('dif8').selection('input').set({'rot3'});
model.component('comp1').geom('geom1').feature('dif8').selection('input2').set({'imp7' 'mov8(2)' 'mov9(2)'});
model.component('comp1').geom('geom1').create('dif4', 'Difference');
model.component('comp1').geom('geom1').feature('dif4').label('Difference 4 (brain - encap)');
model.component('comp1').geom('geom1').feature('dif4').set('keepadd', true);
model.component('comp1').geom('geom1').feature('dif4').set('keepsubtract', true);
model.component('comp1').geom('geom1').feature('dif4').selection('input').set({'imp7'});
model.component('comp1').geom('geom1').feature('dif4').selection('input2').set({'mov8(2)' 'mov9(2)'});
model.component('comp1').geom('geom1').create('dif6', 'Difference');
model.component('comp1').geom('geom1').feature('dif6').label('Difference 6 (nucl - encap)');
model.component('comp1').geom('geom1').feature('dif6').selection('input').set({'imp8(1)' 'imp9(1)' 'imp10(1)'});
model.component('comp1').geom('geom1').feature('dif6').selection('input2').set({'mov8(2)' 'mov9(2)'});
model.component('comp1').geom('geom1').create('del1', 'Delete');
model.component('comp1').geom('geom1').feature('del1').selection('input').init(3);
model.component('comp1').geom('geom1').feature('del1').selection('input').set('sph1(1)', 1);
model.component('comp1').geom('geom1').feature('del1').selection('input').set('rot3(1)', 1);
model.component('comp1').geom('geom1').feature('del1').selection('input').set('imp7(1)', 1);
model.component('comp1').geom('geom1').create('ige1', 'IgnoreEdges');
model.component('comp1').geom('geom1').feature('ige1').label('Ignore Edges');
model.component('comp1').geom('geom1').feature('ige1').selection('input').set('fin(1)', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 64 66 76 78 80 82 84 85 92 96 98 100 110 112 114 116 118 119 126 130 132 134 144 146 148 150 152 153 160 164 166 168 178 180 182 184 186 187 192 193 194 195 196 200 201 202 203 205 206 207 208 209 210 211 212 213 214 215 216 221 226 227 229 231 238 239 240 243 248 249 252 254 264 269 271 273 278 279 280 284 286 290 291 292 300 304 305 306 308 310 311 312 314 319 320 321 325 327 331 335 336 337 341 343 344 345 347 349 352 353 355 356 358 361 362 365 367 369 371 372 373 374 375 377 378 379 380 381]);
model.component('comp1').geom('geom1').create('parf1', 'PartitionFaces');
model.component('comp1').geom('geom1').feature('parf1').selection('face').set('ige1(1)', 6);
model.component('comp1').geom('geom1').feature('parf1').selection('vertexsegment').set('ige1(1)', [3 109]);
model.component('comp1').geom('geom1').create('cmf1', 'CompositeFaces');
model.component('comp1').geom('geom1').feature('cmf1').selection('input').set('parf1(1)', [4 5 8 9 44]);
model.component('comp1').geom('geom1').run;

model.view('view17').tag('view171');
model.view('view2').tag('view17');
model.view('view18').tag('view181');
model.view('view3').tag('view18');
model.view('view19').tag('view191');
model.view('view4').tag('view19');
model.view('view20').tag('view201');
model.view('view5').tag('view20');
model.view('view21').tag('view211');
model.view('view6').tag('view21');
model.view('view22').tag('view221');
model.view('view7').tag('view22');
model.view('view23').tag('view231');
model.view('view8').tag('view23');
model.view('view24').tag('view241');
model.view('view9').tag('view24');
model.view('view25').tag('view251');
model.view('view10').tag('view25');
model.view('view26').tag('view261');
model.view('view11').tag('view26');
model.view('view12').tag('view27');
model.view('view13').tag('view28');
model.view('view14').tag('view29');
model.view('view15').tag('view30');
model.view('view16').tag('view31');
model.view('view171').tag('view32');
model.view('view181').tag('view33');
model.view('view191').tag('view34');
model.view('view201').tag('view35');
model.view('view211').tag('view37');
model.view('view221').tag('view38');
model.view('view231').tag('view39');
model.view('view241').tag('view40');
model.view('view251').tag('view41');
model.view('view261').tag('view42');
model.component('comp1').view('view1').hideEntities.create('hide1');
model.component('comp1').view('view1').hideEntities('hide1').set([1 3]);
model.view.create('view9', 3);
model.view.create('view10', 3);
model.view.create('view11', 3);
model.view.create('view12', 3);
model.view.create('view13', 3);
model.view.create('view14', 3);
model.view.create('view15', 3);
model.view.create('view16', 3);
model.view.create('view36', 3);
model.view.create('view43', 3);

model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('comp1').physics('ec').create('gnd1', 'Ground', 2);
model.component('comp1').physics('ec').feature('gnd1').selection.set([5 6]);
model.component('comp1').physics('ec').create('cucn2', 'CurrentConservation', 3);
model.component('comp1').physics('ec').feature('cucn2').selection.set([4]);
model.component('comp1').physics('ec').create('cucn3', 'CurrentConservation', 3);
model.component('comp1').physics('ec').feature('cucn3').selection.set([1]);
model.component('comp1').physics('ec').create('ncd1', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd1').selection.set([36 37 40]);
model.component('comp1').physics('ec').create('ncd2', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd2').selection.set([33 34 39]);
model.component('comp1').physics('ec').create('ncd3', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd3').selection.set([31 32 38]);
model.component('comp1').physics('ec').create('ncd4', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd4').selection.set([29 30 35]);
model.component('comp1').physics('ec').create('ncd5', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd5').selection.set([23 24 25]);
model.component('comp1').physics('ec').create('ncd6', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd6').selection.set([20 21 22]);
model.component('comp1').physics('ec').create('ncd7', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd7').selection.set([17 18 19]);
model.component('comp1').physics('ec').create('ncd8', 'NormalCurrentDensity', 2);
model.component('comp1').physics('ec').feature('ncd8').selection.set([14 15 16]);

model.component('comp1').mesh('mesh1').autoMeshSize(2);

model.result.table('evl3').label('Evaluation 3D');
model.result.table('evl3').comments('Interactive 3D values');

model.component('comp1').view('view1').set('transparency', true);
model.component('comp1').view('view1').set('uniformblending', true);
model.view('view17').label('View 17');
model.view('view18').label('View 18');
model.view('view19').label('View 19');
model.view('view20').label('View 20');
model.view('view21').label('View 21');
model.view('view22').label('View 22');
model.view('view23').label('View 23');
model.view('view24').label('View 24');
model.view('view25').label('View 25');
model.view('view26').label('View 26');
model.view('view27').label('View 27');
model.view('view28').label('View 28');
model.view('view29').label('View 29');
model.view('view30').label('View 30');
model.view('view31').label('View 31');
model.view('view32').label('View 32');
model.view('view33').label('View 33');
model.view('view34').label('View 34');
model.view('view35').label('View 35');
model.view('view37').label('View 37');
model.view('view38').label('View 38');
model.view('view39').label('View 39');
model.view('view40').label('View 40');
model.view('view41').label('View 41');
model.view('view42').label('View 42');

model.component('comp1').physics('ec').prop('PortSweepSettings').set('PortParamName', 'PortName');
model.component('comp1').physics('ec').feature('cucn1').set('sigma_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn1').set('sigma', {'xx'; 'xy'; 'xz'; 'xy'; 'yy'; 'yz'; 'xz'; 'yz'; 'zz'});
model.component('comp1').physics('ec').feature('cucn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn1').set('epsilonr', {'ee'; '0'; '0'; '0'; 'ee'; '0'; '0'; '0'; 'ee'});
model.component('comp1').physics('ec').feature('cucn1').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec').feature('cucn1').set('minput_numberdensity', 0);
model.component('comp1').physics('ec').feature('dcont1').set('pairDisconnect', true);
model.component('comp1').physics('ec').feature('dcont1').label('Continuity');
model.component('comp1').physics('ec').feature('cucn2').set('sigma_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn2').set('sigma', [0.065095; 0; 0; 0; 0.065095; 0; 0; 0; 0.065095]);
model.component('comp1').physics('ec').feature('cucn2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn2').set('epsilonr', [29790; 0; 0; 0; 29790; 0; 0; 0; 29790]);
model.component('comp1').physics('ec').feature('cucn2').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec').feature('cucn2').set('minput_numberdensity', 0);
model.component('comp1').physics('ec').feature('cucn2').label('Current Conservation 2_encap');
model.component('comp1').physics('ec').feature('cucn3').set('sigma_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn3').set('sigma', [0.3; 0; 0; 0; 0.3; 0; 0; 0; 0.3]);
model.component('comp1').physics('ec').feature('cucn3').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('ec').feature('cucn3').set('epsilonr', [88; 0; 0; 0; 88; 0; 0; 0; 88]);
model.component('comp1').physics('ec').feature('cucn3').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec').feature('cucn3').set('minput_numberdensity', 0);
model.component('comp1').physics('ec').feature('cucn3').label('Current Conservation 3_head');
model.component('comp1').physics('ec').feature('ncd1').set('nJ', '0.4[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd1').label('Normal Current Density 1 (STN_R Contact 1)');
model.component('comp1').physics('ec').feature('ncd2').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd2').label('Normal Current Density 2 (STN_R Contact 2)');
model.component('comp1').physics('ec').feature('ncd3').set('nJ', '0.4[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd3').label('Normal Current Density 3 (STN_R Contact 3)');
model.component('comp1').physics('ec').feature('ncd4').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd4').label('Normal Current Density 4 (STN_R Contact 4)');
model.component('comp1').physics('ec').feature('ncd5').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd5').label('Normal Current Density 5 (GPe_R Contact 1)');
model.component('comp1').physics('ec').feature('ncd6').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd6').label('Normal Current Density 6 (GPe_R Contact 2)');
model.component('comp1').physics('ec').feature('ncd7').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd7').label('Normal Current Density 7 (GPe_R Contact 3)');
model.component('comp1').physics('ec').feature('ncd8').set('nJ', '0[mA]/(0.3514[mm][mm])');
model.component('comp1').physics('ec').feature('ncd8').label('Normal Current Density 8 (GPe_R Contact 4)');

model.study.create('std1');
model.study('std1').create('freq', 'Frequency');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('int2_ds1', 'Grid3D');
model.result.dataset.create('dset2', 'Solution');
model.result.dataset.create('dset3', 'Solution');
model.result.dataset.create('dset4', 'Solution');
model.result.dataset('int2_ds1').set('data', 'none');
model.result.dataset('dset2').selection.geom('geom1', 3);
model.result.dataset('dset2').selection.set([2 4]);
model.result.dataset('dset3').selection.geom('geom1', 3);
model.result.dataset('dset3').selection.set([5 6 7]);
model.result.dataset('dset4').selection.geom('geom1', 2);
model.result.dataset('dset4').selection.set([14 15 16 17 18 19 20 21 22 23 24 25 29 30 31 32 33 34 35 36 37 38 39 40]);
model.result.create('pg1', 'PlotGroup3D');
model.result.create('pg2', 'PlotGroup3D');
model.result.create('pg3', 'PlotGroup3D');
model.result.create('pg4', 'PlotGroup3D');
model.result.create('pg5', 'PlotGroup3D');
model.result('pg1').create('mslc1', 'Multislice');
model.result('pg2').create('iso1', 'Isosurface');
model.result('pg3').create('plot1', 'Slice');
model.result('pg3').feature('plot1').set('data', 'dset1');
model.result('pg3').feature('plot1').set('expr', 'ec.sigmaxx');
model.result('pg4').create('surf1', 'Surface');
model.result('pg4').create('surf2', 'Surface');
model.result('pg4').create('surf3', 'Surface');
model.result('pg4').create('vol1', 'Volume');
model.result('pg4').create('iso1', 'Isosurface');
model.result('pg4').feature('surf1').set('data', 'dset3');
model.result('pg4').feature('surf1').create('tran1', 'Transparency');
model.result('pg4').feature('surf2').set('data', 'dset2');
model.result('pg4').feature('surf2').create('tran1', 'Transparency');
model.result('pg4').feature('surf3').set('data', 'dset4');
model.result('pg4').feature('vol1').set('data', 'dset3');
model.result('pg4').feature('iso1').set('data', 'dset1');
model.result('pg4').feature('iso1').create('tran1', 'Transparency');
model.result('pg5').create('surf1', 'Surface');
model.result('pg5').create('surf2', 'Surface');
model.result('pg5').create('surf3', 'Surface');
model.result('pg5').create('vol1', 'Volume');
model.result('pg5').feature('surf1').set('data', 'dset3');
model.result('pg5').feature('surf1').create('tran1', 'Transparency');
model.result('pg5').feature('surf2').set('data', 'dset2');
model.result('pg5').feature('surf2').create('tran1', 'Transparency');
model.result('pg5').feature('surf3').set('data', 'dset4');
model.result('pg5').feature('vol1').set('data', 'dset3');

model.study('std1').feature('freq').set('plist', 3049);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Frequency Domain');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'3049[Hz]'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').active(true);
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('dDef').set('thresh', 0.1);
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', [3049]);
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('s1').feature('i1').label('Iterative 1.1');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'bicgstab');
model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').set('thresh', 0.1);
model.sol('sol1').runAll;

model.result.dataset('int2_ds1').set('function', 'int2');
model.result.dataset('int2_ds1').set('parmin1', -9.017);
model.result.dataset('int2_ds1').set('parmax1', 80.983);
model.result.dataset('int2_ds1').set('parmin2', -60.4694);
model.result.dataset('int2_ds1').set('parmax2', 53.8786);
model.result.dataset('int2_ds1').set('parmin3', -80.9228);
model.result.dataset('int2_ds1').set('parmax3', 5.3197);
model.result.dataset('dset1').set('frametype', 'material');
model.result.dataset('dset2').label('DBS leads');
model.result.dataset('dset2').set('frametype', 'material');
model.result.dataset('dset3').label('Nuclei');
model.result.dataset('dset3').set('frametype', 'material');
model.result.dataset('dset4').label('Contacts');
model.result.dataset('dset4').set('frametype', 'material');
model.result('pg1').label('Electric Potential (ec)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('mslc1').set('colortable', 'RainbowLightClassic');
model.result('pg1').feature('mslc1').set('resolution', 'normal');
model.result('pg2').set('showlegendsmaxmin', true);
model.result('pg2').set('showlegendsunit', true);
model.result('pg2').feature('iso1').set('levelrounding', false);
model.result('pg2').feature('iso1').set('colortable', 'RainbowClassic');
model.result('pg2').feature('iso1').set('smooth', 'internal');
model.result('pg2').feature('iso1').set('resolution', 'normal');
model.result('pg3').set('titletype', 'manual');
model.result('pg3').set('title', 'ee(x,y,z) (1)');
model.result('pg3').set('edges', false);
model.result('pg3').feature('plot1').set('descr', 'Electrical conductivity, xx component');
model.result('pg3').feature('plot1').set('quickxnumber', 1);
model.result('pg3').feature('plot1').set('interactive', true);
model.result('pg3').feature('plot1').set('shift', -0.00600000000000001);
model.result('pg3').feature('plot1').set('colortable', 'RainbowClassic');
model.result('pg3').feature('plot1').set('smooth', 'none');
model.result('pg3').feature('plot1').set('resolution', 'normal');
model.result('pg4').set('data', 'int2_ds1');
model.result('pg4').set('solrepresentation', 'solnum');
model.result('pg4').set('view', 'view43');
model.result('pg4').set('showhiddenobjects', true);
model.result('pg4').set('edges', false);
model.result('pg4').feature('surf1').set('titletype', 'none');
model.result('pg4').feature('surf1').set('rangecoloractive', true);
model.result('pg4').feature('surf1').set('rangecolormax', 3);
model.result('pg4').feature('surf1').set('colorlegend', false);
model.result('pg4').feature('surf1').set('smooth', 'internal');
model.result('pg4').feature('surf1').set('resolution', 'normal');
model.result('pg4').feature('surf1').feature('tran1').set('uniformblending', 0);
model.result('pg4').feature('surf2').set('smooth', 'internal');
model.result('pg4').feature('surf2').set('allowmaterialsmoothing', false);
model.result('pg4').feature('surf2').set('resolution', 'normal');
model.result('pg4').feature('surf2').feature('tran1').set('transparency', 0.8);
model.result('pg4').feature('surf2').feature('tran1').set('uniformblending', 0);
model.result('pg4').feature('surf3').set('titletype', 'none');
model.result('pg4').feature('surf3').set('colorlegend', false);
model.result('pg4').feature('surf3').set('smooth', 'internal');
model.result('pg4').feature('surf3').set('allowmaterialsmoothing', false);
model.result('pg4').feature('surf3').set('resolution', 'normal');
model.result('pg4').feature('vol1').active(false);
model.result('pg4').feature('vol1').set('smooth', 'internal');
model.result('pg4').feature('vol1').set('resolution', 'normal');
model.result('pg4').feature('iso1').set('descractive', true);
model.result('pg4').feature('iso1').set('descr', 'V');
model.result('pg4').feature('iso1').set('colorlegend', false);
model.result('pg4').feature('iso1').set('smooth', 'internal');
model.result('pg4').feature('iso1').set('resolution', 'normal');
model.result('pg4').feature('iso1').feature('tran1').set('transparency', 0.8);
model.result('pg4').feature('iso1').feature('tran1').set('uniformblending', 0);
model.result('pg5').set('data', 'int2_ds1');
model.result('pg5').set('solrepresentation', 'solnum');
model.result('pg5').set('view', 'view43');
model.result('pg5').set('showhiddenobjects', true);
model.result('pg5').set('edges', false);
model.result('pg5').feature('surf1').set('titletype', 'none');
model.result('pg5').feature('surf1').set('rangecoloractive', true);
model.result('pg5').feature('surf1').set('rangecolormax', 3);
model.result('pg5').feature('surf1').set('colorlegend', false);
model.result('pg5').feature('surf1').set('smooth', 'internal');
model.result('pg5').feature('surf1').set('resolution', 'normal');
model.result('pg5').feature('surf1').feature('tran1').set('uniformblending', 0);
model.result('pg5').feature('surf2').set('smooth', 'internal');
model.result('pg5').feature('surf2').set('allowmaterialsmoothing', false);
model.result('pg5').feature('surf2').set('resolution', 'normal');
model.result('pg5').feature('surf2').feature('tran1').set('transparency', 0.8);
model.result('pg5').feature('surf2').feature('tran1').set('uniformblending', 0);
model.result('pg5').feature('surf3').set('titletype', 'none');
model.result('pg5').feature('surf3').set('colorlegend', false);
model.result('pg5').feature('surf3').set('smooth', 'internal');
model.result('pg5').feature('surf3').set('allowmaterialsmoothing', false);
model.result('pg5').feature('surf3').set('resolution', 'normal');
model.result('pg5').feature('vol1').active(false);
model.result('pg5').feature('vol1').set('smooth', 'internal');
model.result('pg5').feature('vol1').set('resolution', 'normal');

out = model;
