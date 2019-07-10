<?php


echo <<<END

<html>
  
  <head>
    <title>NAM Cf3 Distances Tool</title>
    
    <META http-equiv="Content-Language" content="en-us"/>
    <META http-equiv="Content-Type" content="text/html; charset=UTF-8"/>
    <META name="DC.author" content="Ehsan Kourkchi"/>

    <link rel="stylesheet" href="https://cdn.pydata.org/bokeh/release/bokeh-0.12.16.min.css" type="text/css" />
            
    <script type="text/javascript" src="https://cdn.pydata.org/bokeh/release/bokeh-0.12.16.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>   
    
    
</head>
<body>



END;



if (isset($_POST['function'])) 
	$func=$_POST['function']; 
		else $func="sin(x)";
		
if (isset($_POST['xmin'])) 
	$xmin=$_POST['xmin']; 
		else $xmin='-10';		
		
if (isset($_POST['xmax'])) 
	$xmax=$_POST['xmax']; 
		else $xmax='10';		

$cmd = "bash test.bash ".$func." ".$xmin." ".$xmax;

$command = escapeshellcmd($cmd);
$output = shell_exec($command);

echo '<table>';
// echo '<tr><td>'.$func.'</td></tr>';
echo '<tr><td>'.$output.'</td></tr>';

echo '<tr><td>';
echo '<form method="post" action="ehsan.php">';
echo '  f(x):= ';
echo '  <input type="text" name="function" value='.$func.'>';
echo '  <input type="submit" value="Submit">';
echo '<br> min (x): <input type="text" size="4" maxlength="6" name="xmin" value='.$xmin.' autocomplete="off">';
echo '<br> max (x): <input type="text" size="4" maxlength="6" name="xmax" value='.$xmax.' autocomplete="off">';
echo '</form></td></tr>';

echo '</table>';


echo <<<END

</body>
</html>

END;




