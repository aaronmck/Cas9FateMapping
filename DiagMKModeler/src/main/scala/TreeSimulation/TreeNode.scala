package main.scala.TreeSimulation

/**
 * Our recursive class that walks down the TreeNode structure
 */
case class TreeNode(leftTree: TreeModel, rightTree: TreeModel) {
  def terminalLeft(): Boolean = leftTree.leaves.size < 1
  def terminalRight(): Boolean = rightTree.leaves.size < 1

  var leftBranch :  Option[TreeNode] = None
  var rightBranch : Option[TreeNode] = None

  def branchLeftTree(): Unit = {
    if (!terminalLeft)
      leftBranch = Some(leftTree.splitOnEvent())
    if (!terminalRight)
      rightBranch = Some(leftTree.splitOnEvent())
  }
}

