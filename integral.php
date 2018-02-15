<?php
echo shell_exec(sprintf('/home/gabriel/Desktop/spiroffs-project/calculate-integral %s %s',
  escapeshellarg($_GET['a']),
  escapeshellarg($_GET['b'])
));
?>
