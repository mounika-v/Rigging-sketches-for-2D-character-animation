#include "gl_framework.hpp"
// #include "cg.hpp"

extern GLfloat xrot,yrot,zrot;
extern GLfloat xtran,ytran,ztran;
extern char mode;

namespace csX75
{
  //! Initialize GL State
  void initGL(void)
  {
    //Set framebuffer clear color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    //Set depth buffer furthest depth
    glClearDepth(1.0);
    //Set depth test to less-than
    glDepthFunc(GL_LESS);
    //Enable depth testing
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
  }


  //!GLFW Error Callback
  void error_callback(int error, const char* description)
  {
    std::cerr<<description<<std::endl;
  }

  //!GLFW framebuffer resize callback
  void framebuffer_size_callback(GLFWwindow* window, int width, int height)
  {
    //!Resize the viewport to fit the window size - draw to entire window
    glViewport(0, 0, width, height);
  }

  //!GLFW keyboard callback
  void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
  {
    //!Close the window if the ESC key was pressed
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
      glfwSetWindowShouldClose(window, GL_TRUE);
    else if (key == GLFW_KEY_E && action == GLFW_PRESS && mode == 'I')
    {
      mode = 'E';
      std::cout<<"Now in editing mode. Press M to merge and A to add new. To exit edit mode press I"<<std::endl;
    }
    else if (key == GLFW_KEY_L && action == GLFW_PRESS && ( mode == 'I' || mode == 'S'))
    {
      mode = 'L';
      std::cout<<"add leaf"<<std::endl;
    }
    else if (key == GLFW_KEY_V && action == GLFW_PRESS && ( mode == 'I' || mode == 'S'))
    {
      mode = 'V';
      std::cout<<"Select the vertex to move"<<std::endl;
    }
    else if (key == GLFW_KEY_I && action == GLFW_PRESS)
    {
        mode = 'I';
        std::cout<<"Idle mode activated. Press E to edit"<<std::endl;
    }
    else if (key == GLFW_KEY_M && action == GLFW_PRESS && ( mode == 'I' || mode == 'S'))
    {
        mode = 'M';
        std::cout<<"Merge mode. Select source and destination nodes"<<std::endl;
    }
    else if (key == GLFW_KEY_A && action == GLFW_PRESS && ( mode == 'I' || mode == 'S'))
    {
        mode = 'A';
        std::cout<<"Add mode. Select source, destination and click a point on edge"<<std::endl;
    }
    else if (key == GLFW_KEY_S && action == GLFW_PRESS)
    {
        mode = 'S';
        std::cout<<"Save mode. click anywhere to save"<<std::endl;
    }
    // else if (key == GLFW_KEY_LEFT && action == GLFW_PRESS)
    //   yrot -= 1.0;
    // else if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS)
    //   yrot += 1.0;
    // else if (key == GLFW_KEY_UP && action == GLFW_PRESS)
    //   xrot += 1.0;
    // else if (key == GLFW_KEY_DOWN && action == GLFW_PRESS)
    //   xrot += 1.0;
    // else if (key == GLFW_KEY_PAGE_UP && action == GLFW_PRESS)
    //   zrot += 1.0;
    // else if (key == GLFW_KEY_PAGE_DOWN && action == GLFW_PRESS)
    //   zrot += 1.0;
  }

};
