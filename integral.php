<?php
echo shell_exec(sprintf('/home/gabriel/Desktop/spiroffs-project/calculate-integral %s %s',
  implode(' ', array_map("escapeshellarg", $_GET['a'])),
  implode(' ', array_map("escapeshellarg", $_GET['b']))
));
?>
