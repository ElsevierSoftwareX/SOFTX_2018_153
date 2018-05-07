<?php
header('Content-Type: text/plain');
$code_path = '/home/gabriel/Desktop/spiroffs-project/src/calculate-integral';
echo shell_exec(sprintf("$code_path %s %s",
  implode(' ', array_map("escapeshellarg", $_GET['a'])),
  implode(' ', array_map("escapeshellarg", $_GET['b']))
));
?>
